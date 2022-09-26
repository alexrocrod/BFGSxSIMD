#pragma once

/// Alexandre Rodrigues, 29/07/2021
/// Improved BFGS interface without Debug for ORE

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <LBFGS.h>
#include <vector>
#include <tuple>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using std::make_tuple;

namespace ore {
namespace analytics {
template<typename mmType>
class BfgsOre
{
private:

    uint64_t ncallsR, ncallsF;
    std::shared_ptr<aadc::AADCFunctions<mmType> > aad_funcs; 
    std::shared_ptr<aadc::AADCWorkSpace<mmType> > ws;
    std::vector<aadc::AADCResult> output_args;

    tuple<string, YieldCurveType, string> key;  // curve key
    AADCCurveData curve; // Curve data form the used key  (instead of v_index)

    std::vector<double> ci; // final target outputs
    std::vector<double> out;    // save outputs to calculate gradient
    
    // static const int avx_count = sizeof(mmType) / sizeof(double);
    const int avx_count = sizeof(mmType) / sizeof(double);
    double delta = 1e-7; // default: 1e-7
    const float partial_step = 2 / float(avx_count);

    // Coefficients for derivative calculation:
    const std::vector<double> coeffs_2 = {-1/2.0, 1/2.0};
    const std::vector<double> coeffs_4 = {1/12.0, -2/3.0, 2/3.0, -1/12.0};
    const std::vector<double> coeffs_6 = {-1/60.0, 3/20.0, -3/4.0, 3/4.0, -3/20.0, 1/60.0};
    const std::vector<double> coeffs_8 = {1/280.0, -4/105.0, 1/5.0, -4/5.0, 4/5.0, -1/5.0, 4/105.0, -1/280.0};
    std::vector<double> coeffs;
    
    // PolyFit order (0 disables polyfit):
    int polyfit_order = avx_count - 1;

    // Dg order (only 2, 4, 6 and 8 are valid):
    int dg_count = avx_count;

public:
    BfgsOre(std::shared_ptr<aadc::AADCFunctions<mmType>> aad_funcsi,
    std::shared_ptr<aadc::AADCWorkSpace<mmType> > wsi,
    std::vector<aadc::AADCResult> f_resi,
    std::vector<double> _ci, map<tuple<string, YieldCurveType, string>, AADCCurveData>& _curves,
    tuple<string, YieldCurveType, string> _key = make_tuple("xois_eur", YieldCurveType::Discount, "EUR") ) 
    :   ncallsR(0), ncallsF(0), aad_funcs(aad_funcsi), ws(wsi), output_args(f_resi), key(_key), ci(_ci)
    {
        out = std::vector<double>(output_args.size(), 0.0); // initialize 
        curve = _curves[key];
        resetNcalls();

        coeffs = (avx_count == 8) ? coeffs_8 : coeffs_4;
    }
    // tradicional operator() for LBFGS++ (used before the start of the main loop)
    double operator()(const VectorXd& x, VectorXd& grad)
    {   
        for (int j = 0; j < x.size(); j++) 
            ws->setVal(curve.curve_args[j], x[j]);    

        aad_funcs->forward(*ws);
        ncallsF++;
        
        double fx = 0.0;
        
        for(size_t i = 0; i < output_args.size(); i++) 
        {
            out[i] = ws->valp(output_args[i])[0];
            fx += (out[i] - ci[i]) * (out[i] - ci[i]); 
        }
        
        grad = VectorXd::Zero(x.size());

        for(size_t i = 0; i < output_args.size(); i++) 
        {
            ws->resetDiff();
            ws->setDiff(output_args[i], 1);

            aad_funcs->reverse(*ws);
            ncallsR++;

            for (int j = 0; j < x.size(); j++)
                grad[j] += ws->diffp(curve.curve_args[j])[0] * 2 * (out[i] - ci[i]); 
        }
        return fx; 
    }
    /// operator that computes the fx and gradient in all points given
    double operator()(VectorXd &x, VectorXd& grad, double &step, const VectorXd &drt, const VectorXd &xp)
    {   
        aadc::mmVector<mmType> xs = aadc::mmVector<mmType>(x.size()); // vector of x mmtypes
        VectorXd steps(avx_count);

        for (int i = 0; i < avx_count; i++)
        {
            steps[i] = step * (i + 1) * partial_step;
            for (int j = 0; j < x.size(); j++)
                xs[j][i] = xp[j] + steps[i] * drt[j];
        }

        VectorXd fxs = getForward(xs);
        
        int minElIndex = std::min_element(fxs.begin(), fxs.end()) - fxs.begin(); // index of minimum fx

        double fx = fxs[minElIndex];

        if (polyfit_order > 0)     // calculates the minimum with polyfit:
        { 
            VectorXd coeffs_poly = polyfit(steps, fxs, polyfit_order);

            std::vector<double> zeros = roots(coeffs_poly);

            double stepmin = 4 * step; // maximum acceptable step
            
            for (size_t i = 0; i < zeros.size(); i++)
            {
                if (zeros[i] > 0.5 * steps[0] && zeros[i] < stepmin) 
                    stepmin = zeros[i];
            }

            VectorXd xmin;
            xmin.noalias() = xp + stepmin * drt;

            double fmin = getForward(xmin);

            if (fmin < fx)
            {
                x.noalias() = xmin;
                step = stepmin;
                return fmin;
            }
        }
        
        for (int64_t i = 0; i < x.size(); i++)
            x[i] = aadc::toDblPtr(xs[i])[minElIndex];
        
        step = steps[minElIndex];   // updates the step with the used step

        return fx;
    }
    /// function that calculates the 1D-derivative 
    double getDg(const VectorXd &drt, const VectorXd &x, const double &fx0, VectorXd &grad)   
    {
        if (dg_count < 2 )
        {
            getGrad(x, grad);
            return grad.dot(drt);
        }
                
        aadc::mmVector<mmType> xs(x.size());  // vector of x mmtypes

        int dg2 = dg_count / 2;

        for (int i = 0; i < dg2; i++)
        {
            for (int j = 0; j < x.size(); j++)
            {
                xs[j][i]          = x[j] - (dg2 - i) * delta * drt[j]; // points before x    
                xs[j][i + dg2]    = x[j] + (dg2 - i) * delta * drt[j]; // points after x
            }
        }
    
        VectorXd fx = getForward(xs);

        double res = 0;

        for (int i = 0; i < dg_count; i++)
        {
            res += fx[i] * coeffs[i];
        }

        return res / delta; 
    }
    // function to compute the gradient when required
    void getGrad(const VectorXd& x, VectorXd& grad)
    {   
        for (int64_t j = 0; j < x.size(); j++) 
            ws->setVal(curve.curve_args[j], x[j]);  

        aad_funcs->forward(*ws);
        ncallsF++;

        for(size_t i = 0; i < output_args.size(); i++)
            out[i] = ws->valp(output_args[i])[0];

        grad = VectorXd::Zero(grad.size());

        for(size_t i = 0; i < output_args.size(); i++)
        {
            ws->resetDiff();
            ws->setDiff(output_args[i], 1);

            aad_funcs->reverse(*ws);
            ncallsR++;

            for (size_t j = 0; j < curve.curve_args.size(); j++) 
                grad[j] += ws->diffp(curve.curve_args[j])[0] * 2 * (out[i] - ci[i]); 
        }
        
    }
    // function to compute the function value at x
    double getForward(const VectorXd &x)
    {   
        for (long int j = 0; j < x.size(); j++) 
            ws->setVal(curve.curve_args[j], x[j]);

        aad_funcs->forward(*ws);
        ncallsF++;

        double res = 0.0;

        for(size_t i = 0; i < output_args.size(); i++) {
            out[i] = ws->valp(output_args[i])[0];
            res += (out[i] - ci[i]) * (out[i] - ci[i]); 
        }
        return res;
    }
    // function to compute the function values at the various xs
    VectorXd getForward(const aadc::mmVector<mmType>& x)
    {   
        for (size_t j = 0; j < x.size(); j++) 
            ws->setVal(curve.curve_args[j],x[j]);

        aad_funcs->forward(*ws);
        ncallsF++;

        VectorXd res = VectorXd::Zero(avx_count);

        for(size_t i = 0; i < output_args.size(); i++) 
        {
            double* outs = ws->valp(output_args[i]);
            for (long int j = 0; j < avx_count; j++)
                res[j] += (outs[j] - ci[i]) * (outs[j] - ci[i]);  
        }
        return res;
    }
    // Fit a polynomial to function values depending on the step
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
    std::vector<double> roots(const VectorXd& coeffs)
    {
        int matsz = coeffs.size() - 1;

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

        std::vector<double> vret(matsz);

        for(int i = 0; i < matsz; i++)
            vret[i] = eig(i).real();

        return vret;
    }
    const uint64_t getNcallsR() 
    {
        return ncallsR;
    }
    const uint64_t getNcallsF() 
    {
        return ncallsF;
    }
    const void resetNcalls() 
    {
        ncallsR = 0;
        ncallsF = 0;
    }
    const void setDelta(double d) 
    {
        delta = d;
    }
    const double getDelta() 
    {
        return delta;
    }
    const void setPolyFitOrder(int o)
    {
        polyfit_order = std::min(o, avx_count - 1);

        if (o < 0)
            polyfit_order = 0; //disable
    }
    const int getPolyFitOrder()
    {
        return polyfit_order;
    }
    const void setDgOrder(int ord)
    {
        dg_count = aadc::min(ord,avx_count);

        switch (dg_count)
        {
            case 2:
                coeffs = coeffs_2;
                break;
            case 4:
                coeffs = coeffs_4;
                break;
            case 6:
                coeffs = coeffs_6;
                break;
            case 8:
                coeffs = coeffs_8;
                break;
            default: // invalid dg_count -> disabling finite differences
                dg_count = 0;
                break;
        }
    }
    const int getDgOrder()
    {
        return dg_count;
    }
    const void setDgCoeffs(std::vector<double> cs)
    {
        if (cs.size() != dg_count)
            return;
        
        coeffs = cs;        
    }
    const std::vector<double> setDgCoeffs()
    {
        return coeffs;
    }
};
}
}
