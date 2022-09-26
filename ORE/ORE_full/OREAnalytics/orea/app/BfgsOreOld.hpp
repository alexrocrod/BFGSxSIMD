#pragma once

/// Improved BFGS interface for ORE, only 1 curve for original LBFGS++ implementation

#include <Eigen/Core>

#include <LBFGS.h>
#include <vector>
#include <tuple>

using Eigen::VectorXd;
using std::make_tuple;

namespace ore {
namespace analytics {
template<typename mmType>
class BfgsOreOld
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

public:
    BfgsOreOld(std::shared_ptr<aadc::AADCFunctions<mmType>> aad_funcsi,
    std::shared_ptr<aadc::AADCWorkSpace<mmType> > wsi,
    std::vector<aadc::AADCResult> f_resi,
    std::vector<double> _ci, map<tuple<string, YieldCurveType, string>, AADCCurveData>& _curves,
    tuple<string, YieldCurveType, string> _key = make_tuple("xois_eur", YieldCurveType::Discount, "EUR") ) 
    :   ncallsR(0), ncallsF(0), aad_funcs(aad_funcsi), ws(wsi), output_args(f_resi), key(_key), ci(_ci)
    {
        out = std::vector<double>(output_args.size(), 0.0); // initialize 
        curve = _curves[key];
        resetNcalls();
    }
    double operator()(const VectorXd& x, VectorXd& grad)
    {   
        for (long int j = 0; j < x.size(); j++) 
            ws->setVal(curve.curve_args[j], x[j]);   

        aad_funcs->forward(*ws); 
        ncallsF++;

        double fx = 0.0;

        for(size_t i = 0; i < output_args.size(); i++) 
        {
            out[i] = ws->valp(output_args[i])[0];
            fx += (out[i] - ci[i]) * (out[i] - ci[i]); 
        }

        grad = VectorXd::Zero(grad.size());

        for(size_t i = 0; i < output_args.size(); i++) 
        {
            ws->resetDiff();
            ws->setDiff(output_args[i], 1);

            aad_funcs->reverse(*ws); 
            ncallsR++;

            for (long int j = 0; j < x.size(); j++) 
                grad[j] += ws->diffp(curve.curve_args[j])[0] * 2 * (out[i] - ci[i]); 
        }
        return fx; 
    }
    /// New Operator (Does not use SIMD AVX)
    double operator()(VectorXd &x , VectorXd& grad, double &step, const VectorXd &drt, const VectorXd &xp)
    {   
        return operator()(x,grad); // Tradtional operator
    }
    /// Calculates the directional derivative
    double getDg(const VectorXd &drt, const VectorXd &x, const double &fx0, VectorXd& grad)   
    {
        return grad.dot(drt);
    }
    // function to compute the gradient when required
    void getGrad(const VectorXd& x, VectorXd& grad)
    {   
        return; // grad already calculated in operator()
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
        ncallsF=0;
        ncallsR=0;
    }
};
}
}
