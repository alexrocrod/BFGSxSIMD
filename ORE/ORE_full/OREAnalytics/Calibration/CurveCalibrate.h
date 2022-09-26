#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <LBFGS_old.h>


#include <aadc/aadc_matrix.h>

#include <ql/instruments/makeois.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/time/schedule.hpp>

#include <ql/time/daycounter.hpp>
#include <ql/time/daycounters/all.hpp>
//typedef __m256d mmType;
//#define LENGTH(a) (sizeof(a)/sizeof(a[0]))
#define LENGTH(a) (a.size())
using Eigen::VectorXd;
using namespace LBFGSpp;


namespace aadc {
    template<class T>
    struct VectorType {
        typedef std::vector<T> VecType;
    };

    template<>
    struct VectorType<__m256d> {
        typedef mmVector<__m256d> VecType;
    };

    template<>
    struct VectorType<__m512d> {
        typedef mmVector<__m512d> VecType;
    };
    template<>
    struct VectorType<AADCArgument> {
        typedef VectorArg VecType;
    };

    template<>
    struct VectorType<AADCResult> {
        typedef VectorRes VecType;
    };
}; 


template<typename mmType>
class ScalarSetWSInputsFromVectorVisitor {
public:
    ScalarSetWSInputsFromVectorVisitor(aadc::AADCWorkSpace<mmType>& _ws) : ws(_ws) {}
public:
    void visit(const std::vector<aadc::AADCArgument>& arg_vec, const std::vector<double>& vals) {
        for (int i=0; i<arg_vec.size(); i++) {
            ws.val(arg_vec[i]) = aadc::mmSetConst<mmType>(vals[i]);
        }
    }
private:
    aadc::AADCWorkSpace<mmType>& ws;
};

class InitializeAndMarkAsInputVisitor {
public:
    void visit(const std::vector<idouble>& var_vec, std::vector<aadc::AADCArgument>& arg_vec) const {
        for (int i=0; i<var_vec.size(); i++) {
            arg_vec.push_back(var_vec[i].markAsInput());
        }
    }
};

class InitializeAndCreateDoubleCopyVisitor {
public:
    void visit(const std::vector<idouble>& var_vec, std::vector<double>& dbl_vec) const {
        for (int i=0; i<var_vec.size(); i++) {
            dbl_vec.push_back(var_vec[i].val);
        }
    }
};

template<class VecType>
class InputToFlatVectorVisitor {
typedef typename VecType::value_type ValueType; 
public:
    InputToFlatVectorVisitor (VecType& x_) : x(x_) {}
    template<typename T> void visit(const ValueType& val, const T& val_not_used) {x[i++]=val;}

    
    void visit(const aadc::VectorRes& vec, const aadc::VectorRes& not_used) {
        for (int j=0; j<vec.size(); j++) {
            x[i]=vec[j]; 
            i++;
        }    
    }

public:
    VecType& x;
private:
    int i=0;
};

template<class VecType>
class InputToVectorFromFlatVisitor {
//typedef typename VecType::value_type ValueType; 
public:
    InputToVectorFromFlatVisitor (VecType& x_) : x(x_) {}
    template<typename T> void visit(double& val, const T& val_not_used) {val=x[i++];}

    /*
    void visit(const aadc::VectorRes& vec, const aadc::VectorRes& not_used) {
        for (int j=0; j<vec.size(); j++) {
            x[i]=vec[j]; 
            i++;
        }    
    }
*/
public:
    VecType& x;
private:
    int i=0;
};

class InputToFlatCounter {
public:
    template<typename T1, typename T2>
    void visit(const T1& val_not_used0, const T2& val_not_used) {i++;}

    template<typename T2>
    void visit(const aadc::VectorRes& vec, const T2& val_not_used) {i+=vec.size();}

    int size() {return i;}
private:
    int i=0;
};


using namespace QuantLib;
struct CalibrationSet {
    Integer settlementDays;
    Integer maturity;
    TimeUnit unit;
    Rate marketQuote;
};

template<class numberType>
class DualOISCurveModelOutput {
public:
  
    template<class Visitor, class T2>
    void visit(const Visitor& v,  DualOISCurveModelOutput<T2>& other) {
        v.visit(oisSwaps, other.oisSwaps);
        v.visit(fras, other.fras);
        v.visit(euriborSwaps, other.euriborSwaps);
    }

    template<class Visitor, class T2>
    void visitM(Visitor& v,  DualOISCurveModelOutput<T2>& other) {
        v.visit(oisSwaps, other.oisSwaps);
        v.visit(fras, other.fras);
        v.visit(euriborSwaps, other.euriborSwaps);
    }
public:
    typedef typename aadc::VectorType<numberType>::VecType VecType; 
    VecType oisSwaps;
    VecType euriborSwaps;
    VecType fras;
};



template<class numberType>
class DualOISCurveModelInput {
public:
    template<class Visitor, class Other>
    void visit(Visitor& v, const Other& other) {
        for (int i=0; i<discount_rates.size(); i++) {
            v.visit(discount_rates[i], other);
        }
        for (int i=0; i<forecast_rates.size(); i++) {
            v.visit(forecast_rates[i], other);
        }
    }

public:
    typedef typename aadc::VectorType<numberType>::VecType VecType; 
    VecType discount_rates;
    VecType forecast_rates;
};



class CurveCalibration {
public:
    void calibrateCurve();
    void calibrateCurveNewton();
    virtual void CreateCalibCurve()=0;
    virtual void MarkCurveInputs()=0;
    virtual void CreatePortfolio()=0;
    virtual void RevalueCalibPortfolio()=0;
    virtual void PricingAfterCalibration(int n)=0;
    virtual void InitParameters(VectorXd& x, std::vector<aadc::AADCArgument>& x_args)=0;
    virtual void InitResults(VectorXd& y_vals, std::vector<aadc::AADCResult>& y_args)=0;


    void recordPortfolioRevalue() {
        CreateCalibCurve();
        CreatePortfolio();
        
        aad_funcs.startRecording();
            //{AADDebugStart();}
            MarkCurveInputs();
            RevalueCalibPortfolio();
        aad_funcs.stopRecording();
        //AADDebugStop();

    	std::cout << "Number active to passive conversions: " << CAAD_iVarNumExtractPassive() << std::endl;
	    for (int i = 0; i < CAAD_iVarNumExtractPassive(); ++i) {
    		std::cout << CAAD_iVarGetExtractPassiveLocation(i) << ":" << CAAD_iVarGetExtractPassiveLocationLine(i) << std::endl;
    	}
    }

public:
    aadc::AADCFunctions<mmType> aad_funcs;
    aadc::AADCResult res_arg;

    DualOISCurveModelOutput<Real> portfolio_prices;
    DualOISCurveModelOutput<aadc::AADCArgument> portfolio_arg_prices;
    DualOISCurveModelOutput<double> portfolio_init_prices;

    DualOISCurveModelInput<double> curve_calibrated_inputs;
    std::vector<aadc::AADCArgument> calibrated_x_args;
    VectorXd calibrated_x_vals;
    Eigen::VectorXi indices_cutted_P, indices_cutted_Q;
    int rank;    
    
    VectorXd calibrated_rates;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> calibrated_jacobian, calibrated_jacobian_inverse;
};

class DualOISCurveCalibration : public CurveCalibration {
public:
    DualOISCurveCalibration(Date _valuationDate);
    void CreateCalibCurve();
    void CreatePortfolio();
    void MarkCurveInputs();
    void RevalueCalibPortfolio();
    void PricingAfterCalibration(int n);
    void InitParameters(VectorXd& x, std::vector<aadc::AADCArgument>& x_args);
    void InitResults(VectorXd& y_vals, std::vector<aadc::AADCResult>& y_args);

    Calendar calendar;
    std::vector<double> check_rates_values;
public:    
    Date valuationDate;
    Date settlementDate;
    Integer fixingDays;


public:
    boost::shared_ptr<PiecewiseYieldCurve<Discount, LogLinear>> euribor6MTermStructure;
    boost::shared_ptr<PiecewiseYieldCurve<Discount, LogLinear>> eoniaTermStructure;
    RelinkableHandle<YieldTermStructure> discountingTermStructure;
    RelinkableHandle<YieldTermStructure> forecastingTermStructure;

    
    DualOISCurveModelInput<aadc::AADCArgument> curve_arg_inputs;
    DualOISCurveModelInput<double> curve_init_inputs;

    DualOISCurveModelOutput<aadc::AADCResult> portfolio_arg_res_prices;

    std::vector<CalibrationSet>  deposits, oisSwaps, euriborSwaps, fras;

    std::vector<ext::shared_ptr<OvernightIndexedSwap>> swaps; 
    std::vector<ForwardRateAgreement> fra_s;
    std::vector<VanillaSwap> euribor_swaps;
};
