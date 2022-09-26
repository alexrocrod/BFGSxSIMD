/// Alexandre Rodrigues , 21/7/2022 
/// Final General Class for Improved BFGS interaction

#pragma once

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <LBFGS.h>
#include <aadc/aadc.h>
#include <vector>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using namespace aadc;

template<typename mmType>
class BfgsHelper
{
private:
    uint64_t ncallsR, ncallsF;
    std::shared_ptr<aadc::AADCFunctions<mmType> > aad_funcs; 
    std::shared_ptr<aadc::AADCWorkSpace<mmType> > ws;
    std::vector<aadc::AADCArgument> v_index;
    std::vector<aadc::AADCResult> f_res;

    const int avx_count = sizeof(mmType) / sizeof(double);
    int n;
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
    BfgsHelper (std::shared_ptr<aadc::AADCFunctions<mmType>> aad_funcsi,
    std::shared_ptr<aadc::AADCWorkSpace<mmType>> wsi, std::vector<aadc::AADCArgument> v_indexi,
    std::vector<aadc::AADCResult> f_resi ) 
    : ncallsR(0), ncallsF(0), aad_funcs(aad_funcsi), ws(wsi), v_index(v_indexi),  f_res(f_resi), n(v_indexi.size())
    {
        coeffs = (avx_count == 8) ? coeffs_8 : coeffs_4;
    }
    // tradicional operator() for LBFGS++ (used before the start of the main loop)
    double operator()(const VectorXd& x, VectorXd& grad)
    {   
        for (int i = 0; i < n; i++)
            ws->setVal(v_index[i], x[i]);    

        aad_funcs->forward(*ws);
        ncallsF++;
        
        double fx = 0.0;

        for(size_t i = 0; i < f_res.size(); i++)
            fx += ws->valp(f_res[i])[0];

        grad = VectorXd::Zero(n);

        for (size_t i = 0; i < f_res.size(); i++)
        {
            ws->setDiff(f_res[i], 1);
            
            aad_funcs->reverse(*ws);
            ncallsR++;  
            
            for (int j = 0; j < n; j++)
                grad[j] += ws->diffp(v_index[j])[0]; // may be changed 
        }
        return fx;
    }
    /// operator that computes the fx and gradient in all points given
    double operator()(VectorXd &x, VectorXd& grad, double &step, const VectorXd &drt, const VectorXd &xp)
    {   
        mmVector<mmType> xi = mmVector<mmType>(n);      // vector of x mmtypes
        VectorXd steps(avx_count);                      // new steps for each point

        for (int i = 0; i < avx_count; i++)
        {
            steps[i] = step * (i + 1) * partial_step;
            for (int j = 0; j < n; j++)
                xi[j][i] = xp[j] + steps[i] * drt[j];
        }

        VectorXd fxs = getForward(xi);

        int minElIndex = std::min_element(fxs.begin(), fxs.end()) - fxs.begin();    // index of minimum fx

        double fx = fxs[minElIndex];

        if (polyfit_order > 0)    // calculates the minimum with polyfit:
        {    
            VectorXd coeffs_poly = polyfit(steps, fxs, polyfit_order);

            std::vector<double> zeros = roots(coeffs_poly);

            double stepmin = 4 * step;

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

        // if (dg_count < 2)       // not using finite differences
        //     getGrad(x, grad);
        

        step = steps[minElIndex];           // updates the step with the used step

        for (int i = 0; i < n; i++)
            x[i] = toDblPtr(xi[i])[minElIndex];

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
        
        mmVector<mmType> xi(n);  // vector of x mmtypes

        int dg2 = dg_count / 2;

        for (int i = 0; i < dg2; i++)
        {
            for (int j = 0; j < n; j++)
            {
                xi[j][i]       = x[j] - (dg2 - i) * delta * drt[j];    // points before x  
                xi[j][i + dg2] = x[j] + (dg2 - i) * delta * drt[j];    // points after x
            }
        }
    
        VectorXd fx = getForward(xi);

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
        for (int i = 0; i < n; i++)
            ws->setVal(v_index[i], x[i]);    

        aad_funcs->forward(*ws);
        ncallsF++;

        grad = VectorXd::Zero(n);

        for (int i = 0; i < f_res.size(); i++)
        {
            ws->setDiff(f_res[i], 1);
            
            aad_funcs->reverse(*ws);
            ncallsR++;  
            
            for (int j = 0; j < n; j++)
                grad[j] += ws->diffp(v_index[j])[0]; // may be changed 
        }
    }  
    double getForward(const VectorXd& x)
    {   
        for (int i = 0; i < n; ++i)
            ws->setVal(v_index[i], x[i]);    
        

        aad_funcs->forward(*ws); 
        ncallsF++;
        
        double res = 0.0;

        for (int i = 0; i < f_res.size(); i++){
            double outs = ws->valp(f_res[i])[0];
            res += outs;
        }

        return res;
    } 
    VectorXd getForward( const mmVector<mmType>& x)
    {   
        for (int i = 0; i < n; ++i)
            ws->setVal(v_index[i], x[i]);    
        
        aad_funcs->forward(*ws); 
        ncallsF++;
        
        VectorXd res = VectorXd::Zero(avx_count);

        for (int i = 0; i < f_res.size(); i++)
        {
            double* outs = ws->valp(f_res[i]);
            for (int j = 0; j < avx_count; j++)
                res[j] += outs[j];
        }

        return res;
    } 
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
    const void setDelta(double d)
    {
        delta = d;
    }
    const double getDelta()
    {
        return delta;
    }
    const void setPolyFitOrder(int ord)
    {
        polyfit_order = min(ord, avx_count - 1);

        if (polyfit_order < 0)
            polyfit_order = 0; //disable
    }
    const int getPolyFitOrder()
    {
        return polyfit_order;
    }
    const void setDgOrder(int ord)
    {
        dg_count = min(ord,avx_count);

        if (dg_count == 2)
            coeffs = coeffs_2;
        else if (dg_count == 4)
            coeffs = coeffs_4;
        else if (dg_count == 6)
            coeffs = coeffs_6;
        else if (dg_count == 8)
            coeffs = coeffs_8;  
        else // invalid dg_count -> disabling finite differences
            dg_count = 0;
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