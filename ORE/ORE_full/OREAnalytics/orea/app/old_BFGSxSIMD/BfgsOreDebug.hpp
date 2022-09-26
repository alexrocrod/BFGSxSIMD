#pragma once

/// Alexandre Rodrigues, 29/07/2021
/// Improved BFGS interface with Debug for ORE

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <LBFGS.h>
#include <vector>
#include <tuple>

#define LBFGS_AADC_DEBUG 1 // to print both derivative calculations 

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using std::make_tuple;

namespace ore {
namespace analytics {
class BfgsOreDebug
{
private:
#if AADC_512
    typedef __m512d mmType;
    // Coefficients for derivative calculation:
    std::vector<double> coeffs = {1/280.0,-4/105.0,1/5.0,-4/5.0,4/5.0,-1/5.0,4/105.0,-1/280.0};
#else
    typedef __m256d mmType;
    // Coefficients for derivative calculation:
    std::vector<double> coeffs = {1/12.0,-2/3.0,2/3.0,-1/12.0};
#endif 
    uint64_t ncallsR, ncallsF; // calls statistics
    std::shared_ptr<aadc::AADCWorkSpace<mmType> > ws; // AADC WorkSpace
    std::vector<double> ci; // final target outputs
    std::vector<aadc::AADCResult> output_args; // outputes form the AADC function
    AADCCurveData curve; // Curve data form the used key
    std::shared_ptr<aadc::AADCFunctions<mmType> > aad_funcs; // AADC function
    std::vector<double> aadc_npv;    // save outputs to calculate gradient
    tuple<string, YieldCurveType, string> key; // curve key
    static const int avx_count = sizeof(mmType) / sizeof(double); // number of AVX points
    static const int avx_half = avx_count / 2; 
    
    double delta = 1e-14; // default: 1e-14
    const float partial_step = 2 / float(avx_count); // default: 2/float(avx_count)
    int polyfit_order = avx_count - 1;  // default: avx_count-1 (0 disables polyfit)

public:
    BfgsOreDebug( std::shared_ptr<aadc::AADCWorkSpace<mmType> > _ws, std::vector<double> _ci, std::vector<aadc::AADCResult>& _output_args,
        map<tuple<string, YieldCurveType, string>, AADCCurveData>& _curves, std::shared_ptr<aadc::AADCFunctions<mmType> > _aad_funcs,
        tuple<string, YieldCurveType, string> _key = make_tuple("xois_eur", YieldCurveType::Discount, "EUR") ) :
        ws(_ws), ci(_ci), output_args(_output_args), aad_funcs(_aad_funcs), key(_key)
    {
        aadc_npv = std::vector<double>(_output_args.size(), 0.); // initialize 
        curve = _curves[key];
        resetNcalls();
    }
    // tradicional operator() for LBFGS++ (used before the start of the main loop)
    double operator()(const VectorXd& x, VectorXd& grad)
    {   
        std::cout << "Traditonal operator called, x: " << x.transpose() << std::endl;
        for (long int j = 0; j < x.size(); j++) 
            ws->setVal(curve.curve_args[j], x[j]);    

        aad_funcs->forward(*ws); 
        ncallsF++;
        
        double fx = 0.0;
        for(size_t i = 0; i < output_args.size(); i++) {
            aadc_npv[i] = ws->val(output_args[i])[0];
            fx += (aadc_npv[i] - ci[i]) * (aadc_npv[i] - ci[i]); 
        }
        
        grad = VectorXd::Zero(x.size());
        for(size_t i = 0; i < output_args.size(); i++) {
            ws->resetDiff();
            ws->diff(output_args[i]) = aadc::mmSetConst<mmType>(1);
            aad_funcs->reverse(*ws); 
            ncallsR++;
            for (size_t j = 0; j < curve.curve_args.size(); j++)
                grad[j] += ws->diff(curve.curve_args[j])[0] * 2 * (aadc_npv[i] - ci[i]); 
        }
        std::cout << "Forward Calls: " << ncallsF << ", Reverse Calls: " << ncallsR << std::endl; 
        std::cout << "f(x) = " << fx << ", x = " << x.transpose() << std::endl;
        std::cout << "grad = " << grad.transpose() << std::endl;
        return fx; 
    }
    /// operator that computes the fx and gradient in all points given
    double operator()(VectorXd &x, VectorXd& grad, double &step, const VectorXd &drt, const VectorXd &xp)
    {   
        aadc::mmVector<mmType> xs = aadc::mmVector<mmType>(x.size()); // vector of x mmtypes
        VectorXd steps(avx_count);
        for (int i = 0; i < avx_count; i++){
            steps[i] = step * (i + 1) * partial_step;
            for (int j = 0; j < x.size(); j++)
                xs[j][i] = xp[j] + steps[i] * drt[j];
        }

        VectorXd fxs = getForward(xs);
        int minElIndex = std::min_element(fxs.begin(), fxs.end()) - fxs.begin(); // index of minimum fx

        double fx = fxs[minElIndex];
        for (long int i = 0; i < x.size(); i++)
            x[i] = aadc::toDblPtr(xs[i])[minElIndex];

        std::cout << "New operator called, Forward Calls: " << ncallsF << ", Reverse Calls: " << ncallsR << std::endl; 
        std::cout << "Finalized initial approach, x = " << x.transpose() << ", fx= " << fx << std::endl; 

        if (polyfit_order > 0){     // calculates the minimum with polyfit:
            
            std::cout << "Starting Polyfit, order = " << polyfit_order << std::endl; 
        
            VectorXd coeffs_poly = polyfit(steps, fxs, polyfit_order);

            std::cout << "Polyfit values: " << std::endl; 
            for (long int i = 0; i < coeffs_poly.size(); i++)
                std::cout << "p[" << i << "] = " << coeffs_poly[i] << std::endl;

            auto zeros = roots(coeffs_poly);

            double stepmin = 4 * step; // maximum acceptable step
            
            std::cout << "Roots of the polynomial: " << std::endl; 
            for (size_t i = 0; i < zeros.size(); i++){
                std::cout << "roots[" << i << "] = " << zeros[i] << std::endl;
                double zero = zeros[i].real();
                if (zero > 0.5 * steps[0] && zero < stepmin) 
                    stepmin = zero;
            }

            VectorXd xmin;
            xmin.noalias() = xp + stepmin * drt;

            double fmin = getForward(xmin);

            std::cout << "Step choosen with Polyfit: stepmin = " << stepmin  << ", f(xmin) = " << fmin << ", xmin = " << xmin.transpose() << std::endl;

            if (fmin < fx){
                x.noalias() = xmin;
                step = stepmin;
                std::cout << "Used the PolyFit results.\n";
                return fmin;
            }
        }
        step = steps[minElIndex];
        std::cout << "Used the Initial Approach results.\n";
        std::cout << "step = " << step  << ", f(x) = " << fx << ", x = " << x.transpose() << std::endl;


        return fx;
    }
    /// function that calculates the 1D-derivative 
    double getDg(const VectorXd &drt, const VectorXd &x, const double &fx0, VectorXd &grad)   
    {
        std::cout << "getDg() called, delta: " << delta << ", drt: " << drt.transpose() << std::endl; 
        std::cout << "x: " << x.transpose() << std::endl; 
        aadc::mmVector<mmType> xs = aadc::mmVector<mmType>(x.size());  // vector of x mmtypes

        for (int i = 0; i < avx_half; i++)
            for (int j = 0; j < x.size(); j++){
                xs[j][i]          = x[j] - (avx_half-i) * delta * drt[j]; // points before x    
                xs[j][i+avx_half] = x[j] + (avx_half-i) * delta * drt[j]; // points after x
            }
    
        VectorXd fx = getForward(xs);

        double res = 0;
        for (int i = 0; i < avx_count; i++){
            std::cout << "fx[" << i << "] = " << fx[i] << std::endl; 
            res += fx[i] * coeffs[i];
        }
        res /= delta; // directional derivative at x

        std::cout << "Derivative with Finite Diffs: " << res << std::endl;
# if LBFGS_AADC_DEBUG
        getGrad(x, grad); // computes the gradient as it was not done in operator()
        std::cout << "Derivative with Gradient: " << grad.dot(drt) << std::endl;
#endif
        return res; 
    }
    // function to compute the gradient when required
    void getGrad(const VectorXd& x, VectorXd& grad)
    {   
        std::cout << "getGrad() called, x: " << x.transpose() << std::endl; 
        
        for (long int j = 0; j < x.size(); j++) 
            ws->val(curve.curve_args[j]) = aadc::mmSetConst<mmType>(x[j]);

        aad_funcs->forward(*ws); 
        ncallsF++;

        for(size_t i = 0; i < output_args.size(); i++)
            aadc_npv[i] = ws->val(output_args[i])[0];

        grad = VectorXd::Zero(grad.size());
        for(size_t i = 0; i < output_args.size(); i++) {
            ws->resetDiff();
            ws->diff(output_args[i]) = aadc::mmSetConst<mmType>(1);
            aad_funcs->reverse(*ws); 
            ncallsR++;
            for (size_t j = 0; j < curve.curve_args.size(); j++) 
                grad[j] += ws->diff(curve.curve_args[j])[0] * 2 * (aadc_npv[i] - ci[i]); 
        }
        std::cout << "grad: " << grad.transpose() << std::endl; 
        
    }
    // function to compute the function value at x
    double getForward(const VectorXd &x)
    {   
        std::cout << "getForward(VectorXd) called, x: " << x.transpose() << std::endl; 
        for (long int j = 0; j < x.size(); j++) 
            ws->val(curve.curve_args[j]) = aadc::mmSetConst<mmType>(x[j]);

        aad_funcs->forward(*ws); 
        ncallsF++;

        double res = 0.0;

        for(size_t i = 0; i < output_args.size(); i++) {
            aadc_npv[i] = ws->val(output_args[i])[0];
            res += (aadc_npv[i] - ci[i]) * (aadc_npv[i] - ci[i]); 
        }
        std::cout << "f(x) = " << res << std::endl;
        
        return res;
    }
    // function to compute the function values at the various xs
    VectorXd getForward(const aadc::mmVector<mmType>& x)
    {   
        std::cout << "getForward(mmVector) called" << std::endl; 

        for (size_t j=0; j< x.size(); j++) 
            ws->val(curve.curve_args[j]) = x[j];

        aad_funcs->forward(*ws); 
        ncallsF++;

        VectorXd res = VectorXd::Zero(avx_count);

        for(size_t i = 0; i < output_args.size(); i++) {
            double* outs = ws->valp(output_args[i]);
            for (long int j = 0; j < avx_count; j++)
                res[j]+=(outs[j]-ci[i])*(outs[j]-ci[i]);  
        }
        std::cout << "f(xs) = " << res.transpose() << std::endl;

        return res;
    }
    // Fit a polynomial to function values depending on the step
    // from: ...
    VectorXd polyfit(VectorXd xvals, VectorXd yvals, int order) 
    {
        Eigen::MatrixXd A(xvals.size(), order + 1);

        for (int i = 0; i < xvals.size(); i++)
            A(i, 0) = 1.0;

        for (int j = 0; j < xvals.size(); j++) {
            for (int i = 0; i < order; i++)
                A(j, i + 1) = A(j, i) * xvals(j);
        }

        auto Q = A.householderQr();
        auto result = Q.solve(yvals);
        return result;
    }
    // Find the complex roots of the polynomial
    // from: ...
    std::vector<std::complex<double>> roots(const VectorXd& coeffs)
    {
        int matsz = coeffs.size() - 1;
        std::vector<std::complex<double>> vret(matsz);

        MatrixXd companion_mat(matsz,matsz);

        for(int n = 0; n < matsz; n++){
            for(int m = 0; m < matsz; m++){
                if(n == m + 1)
                    companion_mat(n,m) = 1.0;
                if(m == matsz-1)
                    companion_mat(n,m) = -coeffs[n] / coeffs[matsz];
            }
        }

        MatrixXcd eig = companion_mat.eigenvalues();

        for(int i = 0; i < matsz; i++)
            vret[i] = eig(i);

        return vret;
    }
    // Getters and Setters:
    const uint64_t getNcallsR() {
        return ncallsR;
    }
    const uint64_t getNcallsF() {
        return ncallsF;
    }
    const void resetNcalls() {
        ncallsR = 0;
        ncallsF = 0;
    }
    const void setDelta(double d) {
        delta = d;
    }
    const double getDelta() {
        return delta;
    }
    const void setPolyFitOrder(int o){
        polyfit_order = std::min(o, avx_count-1);
        if (o < 0)
            polyfit_order = 0; //disable
    }
    const int getPolyFitOrder(){
        return polyfit_order;
    }
};
}
}
