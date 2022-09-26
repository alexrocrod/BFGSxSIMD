//#include "aadc/aadc_debug.h"

#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <ctime>
#include <random>
#include <aadc/aadc.h>
// #include <aadc/aadc_private.h>
#include <aadc/aadc_matrix.h>



#include <ql/quantlib.hpp>
#include "CurveCalibrate.h"
const int shift_ois=2;
const int shift_eur=1;


namespace QuantLib {
    namespace AAD {
        typedef QuantLib::Calendar Calendar;
        typedef QuantLib::UnitedKingdom UnitedKingdom;
        typedef QuantLib::UnitedStates UnitedStates;
    }
}

using namespace QuantLib;
using namespace Eigen;

DualOISCurveCalibration::DualOISCurveCalibration(Date _valuationDate) : valuationDate(_valuationDate)
{
    calendar = JointCalendar(UnitedKingdom(UnitedKingdom::Exchange), UnitedStates(UnitedStates::Settlement), JoinHolidays);
    Integer fixingDays = 2;
    settlementDate = calendar.advance(valuationDate, fixingDays, Days);
};


void DualOISCurveCalibration::CreateCalibCurve()
{
    calendar = JointCalendar(UnitedKingdom(UnitedKingdom::Exchange), UnitedStates(UnitedStates::Settlement), JoinHolidays);


    /***********************
     ***    CALIBRATION SET    ***
     ***********************/

    deposits ={
        { 0, 1, Days, 1.10 },
        { 1, 1, Days, 1.10 },
        { 2, 1, Days, 1.40 },
    };
    
    oisSwaps ={
//       { 0, 1, Days, 0.01 },
//       { 1, 1, Days, 0.0110 },
       { 2, 1, Weeks, 0.0140 },
       { 2, 2, Weeks, 0.0140 },
       { 2, 3, Weeks, 0.0150 },
       { 2, 1, Months, 0.0170 },
       { 2, 2, Months, 0.0190 },
       { 2, 3, Months, 0.0205 },
       { 2, 4, Months, 0.0208 },
       { 2, 5, Months, 0.0211 },
       { 2, 6, Months, 0.0213 }
    };

    fras =  {
      /* { 0, 1, Months, 0.045 },
         { 0, 2, Months, 0.042720 },
        { 0, 3, Months, 0.042600 },
       { 0, 4, Months, 0.042560 },
        { 0, 5, Months, 0.042520 },
        { 0, 6, Months, 0.042480 },
        { 0, 7, Months, 0.042540 },
        { 0, 8, Months, 0.042610 },
        { 0, 9, Months, 0.042670 },
        { 0, 10, Months, 0.042790 },
        { 0, 11, Months, 0.042910 },
        { 0, 12, Months, 0.043030 },
        { 0, 13, Months, 0.043180 },
        { 0, 14, Months, 0.043350 },
        { 0, 15, Months, 0.043520 },
        { 0, 16, Months, 0.043710 },
        { 0, 17, Months, 0.043890 },
        { 0, 18, Months, 0.044090 }, */ 
    };


  euriborSwaps = {
        { 0, 3, Years,  0.05}, //0.004240 },
        { 0, 4, Years, 0.04760 },
        { 0, 5, Years, 0.047620 },
        { 0, 6, Years, 0.049540 },
        { 0, 7, Years, 0.041350 },
        { 0, 8, Years, 0.043030 },
        { 0, 9, Years, 0.055840 },
        { 0, 10, Years, 0.058090 },
        { 0, 12, Years, 0.050370 },
        { 0, 15, Years, 0.051870 },
        { 0, 20, Years, 0.052340 },
        { 0, 30, Years, 0.052560 },
        { 0, 35, Years, 0.052950 },
        
        
   //     { 0, 40, Years, 0.053480 },
   //     { 0, 45, Years, 0.024210 },
   //     { 0, 50, Years, 0.024630 },
   //     { 0, 60, Years, 0.02463001 },
    };
    /********************
     ***    QUOTES    ***
     ********************/

    Real flatRate = 3;
    ext::shared_ptr<Quote> flatQuote(new SimpleQuote(flatRate)); 

    /*********************
     ***  RATE HELPERS ***
     *********************/
   
    // Deposits
    DayCounter depositDayCounter = Actual360();
    std::vector<ext::shared_ptr<RateHelper> > eoniaInstruments;

    for (Size i = 0; i < LENGTH(deposits); i++) {   
        Period term = deposits[i].maturity * deposits[i].unit;
        ext::shared_ptr<RateHelper> depositHelper(new DepositRateHelper(
            Handle<Quote>(flatQuote),
            term, 3,
            calendar, Following,
            false, depositDayCounter));

        eoniaInstruments.push_back(depositHelper);
    }

    // OIS
    ext::shared_ptr<Eonia> eonia(new Eonia);

    for (Size i = 0; i < LENGTH(oisSwaps); i++) {
        Real rate = 0.01 * oisSwaps[i].marketQuote;      
        Period term = oisSwaps[i].maturity * oisSwaps[i].unit;
        ext::shared_ptr<RateHelper> eoniaHelper(new
            OISRateHelper(oisSwaps[i].settlementDays,
                term,
                Handle<Quote>(flatQuote),
                eonia,
                Handle<YieldTermStructure>(),
                true,
                0));
        eoniaInstruments.push_back(eoniaHelper);
    }


    /*********************
    **  CURVE BUILDING **
    *********************/

    /*********************
    **   EONIA CURVE    **
    *********************/

    DayCounter termStructureDayCounter = Actual365Fixed();
    eoniaTermStructure=boost::make_shared<PiecewiseYieldCurve<Discount, LogLinear>>(
            settlementDate, eoniaInstruments, termStructureDayCounter, 1.0e-15
            );
    eoniaTermStructure->enableExtrapolation();
    discountingTermStructure.linkTo(eoniaTermStructure);


    /*********************
    **    EURIBOR 6M    **
    *********************/
    ext::shared_ptr<IborIndex> euribor6M(new Euribor6M);

    // SETUP SWAPS

    Frequency swFixedLegFrequency = Annual;
    BusinessDayConvention swFixedLegConvention = Unadjusted;
    DayCounter swFixedLegDayCounter = Thirty360(Thirty360::European);
    ext::shared_ptr<IborIndex> swFloatingLegIndex(new Euribor6M);

    // Euribor 6M curve
    std::vector<ext::shared_ptr<RateHelper> > euribor6MInstruments;  
    bool useIndexedFra = true;  

    // Create FRA helpers
    for (Size i = 0; i < LENGTH(fras); i++) {
        Real rate = 0.01 * fras[i].marketQuote;        
        Period term = fras[i].maturity * fras[i].unit;
        ext::shared_ptr<RateHelper> fraHelper(new FraRateHelper(Handle<Quote>(flatQuote),
            fras[i].maturity, euribor6M,
            Pillar::LastRelevantDate, Date()
            //, useIndexedFra));
        ));

        euribor6MInstruments.push_back(fraHelper);

        if (term <= 2 * Days)
        {
            ext::shared_ptr<RateHelper> depositHelper(new DepositRateHelper(
                Handle<Quote>(flatQuote),
                term, 3,
                calendar, Following,
                false, depositDayCounter));
            euribor6MInstruments.push_back(depositHelper);
        }            
    }

    // Create Swap helpers
    for (Size i = 0; i < LENGTH(euriborSwaps); i++) {
        Period maturity = euriborSwaps[i].maturity * 
        euriborSwaps[i].unit;
        ext::shared_ptr<RateHelper> swapHelper(new SwapRateHelper(
            Handle<Quote>(flatQuote), maturity,
            calendar, swFixedLegFrequency,
            swFixedLegConvention, swFixedLegDayCounter,
            swFloatingLegIndex,
            Handle<Quote>(), 0 * Days, discountingTermStructure));

        euribor6MInstruments.push_back(swapHelper);

    }

    euribor6MTermStructure=boost::make_shared<PiecewiseYieldCurve
    <Discount, LogLinear>>(
            settlementDate, euribor6MInstruments, 
            termStructureDayCounter, 1.0e-15
            );

    // Discounting and projection curve rates
   
    eoniaTermStructure->data();
    euribor6MTermStructure->data();

    discountingTermStructure.linkTo(eoniaTermStructure);
    forecastingTermStructure.linkTo(euribor6MTermStructure);
}


void DualOISCurveCalibration::MarkCurveInputs () {

    std::vector<Real>& discountCurveRates(const_cast<std::vector<Real>&>(eoniaTermStructure->data()));
    std::vector<Real>& projectionCurveRates(const_cast<std::vector<Real>&>(euribor6MTermStructure->data()));

    std::cout << "DiscCurveSise" << discountCurveRates.size() <<"\n";
	for (int i=1; i<discountCurveRates.size(); i++) {
        curve_arg_inputs.discount_rates.push_back(discountCurveRates[i].markAsInput());
        curve_init_inputs.discount_rates.push_back(discountCurveRates[i].val);
    }
    std::cout << "PrCurveSise" << projectionCurveRates.size() <<"\n";
	for (int i=1; i<projectionCurveRates.size(); i++) {
        curve_arg_inputs.forecast_rates.push_back(projectionCurveRates[i].markAsInput());
        curve_init_inputs.forecast_rates.push_back(projectionCurveRates[i].val);
    }

    std::cout <<"\n";
	//eoniaTermStructure->interpolationUpdate();
	//euribor6MTermStructure->interpolationUpdate();
    for (Size i = 0; i < LENGTH(oisSwaps)-shift_ois; i++) {
        portfolio_prices.oisSwaps.push_back(oisSwaps[i].marketQuote);  
    }
    for (Size i = 0; i < fras.size(); i++) {
        portfolio_prices.fras.push_back(fras[i].marketQuote); 
    }
    for (Size i = 0; i < LENGTH(euriborSwaps)-shift_eur; i++) {
        portfolio_prices.euriborSwaps.push_back(euriborSwaps[i].marketQuote); 
    }        
    
    portfolio_prices.visit(InitializeAndCreateDoubleCopyVisitor(), portfolio_init_prices);
    portfolio_prices.visit(InitializeAndMarkAsInputVisitor(), portfolio_arg_prices);
}

void DualOISCurveCalibration::CreatePortfolio() {
    bool useIndexedFra = true;  
    ext::shared_ptr<Eonia> eoniaIndex(new Eonia(discountingTermStructure));

    for (Size i = 0; i < LENGTH(oisSwaps)-shift_ois; i++) {
        Period term = oisSwaps[i].maturity * oisSwaps[i].unit;
        Date effectiveDate = settlementDate + oisSwaps[i].settlementDays;
        ext::shared_ptr<OvernightIndexedSwap> swap = MakeOIS(term, eoniaIndex, oisSwaps[i].marketQuote)
            .withEffectiveDate(Null<Date>())
            .withOvernightLegSpread(0.0)
            .withNominal(1.0)
            .withPaymentLag(0.0)
            .withDiscountingTermStructure(discountingTermStructure)
            .withTelescopicValueDates(false);

        swaps.push_back(swap);
    };


    ext::shared_ptr<IborIndex> euriborIndex(new Euribor6M(forecastingTermStructure));
    for (Size i = 0; i < fras.size(); i++) {
        Real rate = 0.01 * fras[i].marketQuote;
        Date maturity = settlementDate + fras[i].maturity * fras[i].unit;   
        std::cout<< "maturity " << maturity << "\n";    
        ForwardRateAgreement fra(settlementDate,
            maturity,
            Position::Long,
            0,
            1,
            euriborIndex,
            discountingTermStructure
        );
          //  ,useIndexedFra);
        fra_s.push_back(fra);
    };

    Frequency fixedLegFrequency = Annual;
    BusinessDayConvention fixedLegConvention = Unadjusted;
    BusinessDayConvention floatingLegConvention = ModifiedFollowing;
    DayCounter fixedLegDayCounter = Thirty360(Thirty360::European);
    Rate fixedRate = 0.003;
    DayCounter floatingLegDayCounter = Actual360();

    // floating leg
    Frequency floatingLegFrequency = Semiannual;
    Real nominal = 1.0;
    VanillaSwap::Type swapType = VanillaSwap::Payer;

    //pricingEngine=boost::make_shared<PricingEngine>(
      //  new DiscountingSwapEngine(discountingTermStructure));
    
    for (Size i = 0; i < LENGTH(euriborSwaps)-shift_eur; i++) {
        Date maturity = settlementDate + euriborSwaps[i].maturity * euriborSwaps[i].unit;
        Spread spread = 0.0;
        Schedule fixedSchedule(settlementDate, maturity,
            Period(fixedLegFrequency),
            calendar, fixedLegConvention,
            fixedLegConvention,
            DateGeneration::Forward, false);

        Schedule floatSchedule(settlementDate, maturity,
            Period(floatingLegFrequency),
            calendar, floatingLegConvention,
            floatingLegConvention,
            DateGeneration::Forward, false);

        VanillaSwap swap(swapType, nominal,
            fixedSchedule, fixedRate, fixedLegDayCounter,
            floatSchedule, euriborIndex, spread,
            floatingLegDayCounter);

        euribor_swaps.push_back(swap);
    };
}


void DualOISCurveCalibration::RevalueCalibPortfolio () {
    bool useIndexedFra = true;  
    Real total_error=0;
    //total_error=AADC_PRINT(total_error);



    for (Size i = 0; i < swaps.size(); i++) {
        //Real NPV = swap->NPV();
        Rate fairRate = swaps[i]->fairRate();
        Rate error = std::sqrt(1)*(portfolio_prices.oisSwaps[i] - fairRate);        
        std::cout << "OIS swap calibration target value: " << error.val << std::endl;
        portfolio_arg_res_prices.oisSwaps.push_back(error.markAsOutput());
        total_error+=1*(portfolio_prices.oisSwaps[i] - fairRate)*(portfolio_prices.oisSwaps[i] - fairRate);
    };
   
    // Run calibration targets for fras   
    for (Size i = 0; i < fra_s.size(); i++) {
        Rate fairRate = fra_s[i].forwardRate();
        Rate error = std::sqrt(0.1)*(portfolio_prices.fras[i]-fairRate );       
        std::cout << "FRA calibration target value: " << error.val << " " << fairRate.val << std::endl;
        portfolio_arg_res_prices.fras.push_back(error.markAsOutput());
        total_error+=1*(fairRate - portfolio_prices.fras[i])*(fairRate - portfolio_prices.fras[i]);
    }    
    
    boost::shared_ptr<PricingEngine>  pricingEngine(new DiscountingSwapEngine(discountingTermStructure));
    

    Rate fairRate;
    for (Size i = 0; i < euribor_swaps.size(); i++) {
        euribor_swaps[i].setPricingEngine(pricingEngine);
        //Real NPV = swap.NPV();
        fairRate = euribor_swaps[i].fairRate();
        //fairRate=AADC_PRINT(fairRate);
        Rate error = ( portfolio_prices.euriborSwaps[i]-fairRate );
        total_error+=(fairRate - portfolio_prices.euriborSwaps[i])*(fairRate - portfolio_prices.euriborSwaps[i]);
        std::cout << "Swap calibration target value: " << error.val << std::endl;
        std::cout << "Fair rate: " << fairRate.val << std::endl;

        portfolio_arg_res_prices.euriborSwaps.push_back(error.markAsOutput());
        //portfolio_idbl_prices.euriborSwaps.push_back(fairRate);
        //portfolio_arg_res_prices.euriborSwaps.push_back(portfolio_idbl_prices.euriborSwaps.back().markAsOutput());
    };
    //total_error-=0.9*(fairRate - portfolio_prices.euriborSwaps.back())*(fairRate - portfolio_prices.euriborSwaps.back());
    std::cout << "Number of instruments: " << LENGTH(euriborSwaps) + LENGTH(fras) + LENGTH(oisSwaps) << "\n";
    //total_error=AADC_PRINT(total_error);
       
    std::cout << "Total error recorded: " << total_error.val << std::endl;
    res_arg=total_error.markAsOutput();
};


void DualOISCurveCalibration::InitParameters (VectorXd& x, std::vector<aadc::AADCArgument>& x_args) {   
    InputToFlatCounter flat_counter; 
    curve_arg_inputs.visit(flat_counter, false);
    x.resize(flat_counter.size());
    x_args.resize(flat_counter.size());
    
    InputToFlatVectorVisitor<VectorXd> to_flat_val_visitor(x);
    InputToFlatVectorVisitor<std::vector<aadc::AADCArgument>> to_flat_arg_visitor(x_args);

    curve_arg_inputs.visit(to_flat_arg_visitor, false);
    curve_init_inputs.visit(to_flat_val_visitor, false);

    curve_calibrated_inputs=curve_init_inputs;

}

void DualOISCurveCalibration::InitResults (VectorXd& y_vals, std::vector<aadc::AADCResult>& y_args) {   
    InputToFlatCounter flat_counter1;
    portfolio_arg_res_prices.visitM(flat_counter1, portfolio_arg_res_prices);
    y_args.resize(flat_counter1.size());
    y_vals.resize(y_args.size());
    InputToFlatVectorVisitor<std::vector<aadc::AADCResult>> to_flat_arg_visitor(y_args);
    portfolio_arg_res_prices.visitM(to_flat_arg_visitor, portfolio_arg_res_prices);
}

void CurveCalibration::calibrateCurveNewton()
{
    // optimisation code
    std::cout << "Start optmisation: " << std::endl;
    std::shared_ptr<aadc::AADCWorkSpace<mmType>> ws(aad_funcs.createWorkSpace());
    VectorXd x_vals ,y_vals;
    std::vector<aadc::AADCArgument> x_args;
    std::vector<aadc::AADCResult> y_args;
    
    InitParameters(x_vals, x_args);
    InitResults(y_vals, y_args);

    Eigen::Matrix<double, Dynamic, Dynamic> jacobian;
    jacobian.resize(y_args.size(), x_args.size());
    
    std::cout << "y " << y_args.size() << " x " << x_args.size() <<"\n";
    
    auto newton_func = [&] (const VectorXd& x,  VectorXd& y_vals_, Eigen::Matrix<double, Dynamic, Dynamic>&  jacobian_)  {
		for (int i=0; i<x.size(); i++) {
			ws->val(x_args[i])=aadc::mmSetConst<mmType>(x[i]);
		}
        ScalarSetWSInputsFromVectorVisitor<mmType> visitor(*ws);
        portfolio_arg_prices.visitM(visitor, portfolio_init_prices);
        //AADDebugStart();
		aad_funcs.forward(*ws);
        std::cout << ws->val(res_arg)[0] << " val\n";
		double t_error=0;
        for (int i=0; i<y_vals.size(); i++) {
			y_vals_[i]=ws->valp(y_args[i])[0];
            //std::cout << i << " " << y_vals_[i]<<"\n";
            t_error+=	y_vals_[i]*	y_vals_[i];
		}
        std::cout << "total square Error " << t_error <<"\n";

        for (int j=0; j<y_args.size(); j++) {
            ws->resetDiff(); 
            ws->diff(y_args[j])=aadc::mmSetConst<mmType>(1.);
            aad_funcs.reverse(*ws);        
            //AADDebugStop();
            for (int i=0; i<x.size(); i++) {
                jacobian_(j,i)=ws->diffp(x_args[i])[0];
                //if (j%20==0) std::cout << j << " diff " << ws->diff(x_args[i])[0] << "\n";
            }
        }
        Eigen::FullPivLU<Eigen::MatrixXd> l(jacobian);
        std::cout << l.rank() << " matrix rank\n";

	};
    
    //Start of calibration   
    std::cout << "x = \n" << x_vals.transpose() << std::endl;
    for (int i=0; i<1; i++) {
        newton_func(x_vals, y_vals, jacobian);
        VectorXd delta_x = jacobian.colPivHouseholderQr().solve(y_vals);
        x_vals=x_vals-delta_x;
    }
    //End of calibration
    std::cout << "x = \n" << x_vals.transpose() << std::endl;
    calibrated_rates.resize(x_vals.size());
    for (int i=0; i<x_vals.size(); i++) {
        calibrated_x_args.push_back(x_args[i]);
        calibrated_rates[i]=x_vals[i];
        std::cout << x_vals[i] <<", "; 
    }


    std::cout <<"\n From the flat x-> structure << curve_calibrated_inputs.forecast_rates.size() " << " \n";

    InputToVectorFromFlatVisitor<VectorXd> from_flat_visitor(x_vals);
    curve_calibrated_inputs.visit(from_flat_visitor, false);

    for(int i=0; i<curve_calibrated_inputs.discount_rates.size(); i++)
        std::cout << curve_calibrated_inputs.discount_rates[i] << " ";
    std::cout <<"\n";
    for(int i=0; i<curve_calibrated_inputs.forecast_rates.size(); i++)
        std::cout << curve_calibrated_inputs.forecast_rates[i] << " ";


    //Record of data to CalibCurve class
    //calibrated_rates=std::make_shared<VectorXd>();
    std::cout << "\n" << x_vals.size() << " " << x_vals.rows() <<" " << jacobian.rows() <<"end \n";
    calibrated_rates.resize(x_vals.size());
    calibrated_rates=x_vals;
    calibrated_jacobian.resize(jacobian.rows(),jacobian.cols());
    calibrated_jacobian=jacobian;
    Eigen::FullPivLU<Eigen::MatrixXd> l(calibrated_jacobian);
    std::cout << l.isInvertible() << " initial matrix is invertible?\n";
    std::cout << l.rank() << " init. matrix rank\n";
    std::cout << l.threshold() << " thres\n";
    rank=l.rank();

    indices_cutted_P.resize(l.rank());
    indices_cutted_Q.resize(l.rank());
    Eigen::VectorXi P_indices=l.permutationP().indices();
    Eigen::VectorXi Q_indices=l.permutationQ().indices();
    
    //std::cout << "\nPErm " << l.permutationP().toDenseMatrix() << "\n";

    for (int i=0; i<l.rank(); i++) {
            indices_cutted_P[i]=P_indices[i];
            indices_cutted_Q[i]=Q_indices[i];
    }

    std::cout << "\n IndP= " <<P_indices.size() << "; IndQ= " << Q_indices.size() << "\n";

    
Eigen::Matrix<double, 5, 3> m;
m << 0, 0, 0,
     0, 0, 0,
     0, 1, 2,
     0, 1, 1,
     0, 0, 0;
 Eigen::FullPivLU<Eigen::MatrixXd> ll(m);
    Eigen::VectorXi P_i=ll.permutationP().indices();
    Eigen::VectorXi Q_i=ll.permutationQ().indices();
std::cout <<"aaaa\n" << P_i <<"\n";
std::cout <<"aaaa\n" << Q_i <<"\n";
std::cout <<"MMMMMMM\n" <<"\n";
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> 
        MM(m({2,3}, {2,1}));
    Eigen::FullPivLU<Eigen::MatrixXd> ll2(MM);
    std::cout << ll2.rank() << " rank " << ll2.rows() << " " << ll2.cols() << "\n";



    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> 
        M(calibrated_jacobian(indices_cutted_P, indices_cutted_Q));
    Eigen::FullPivLU<Eigen::MatrixXd> l2(M);
    std::cout << l2.rank() << " rank " << l2.rows() << " " << l2.cols() << "\n";
    calibrated_jacobian_inverse.resize(l2.rank(), l2.rank());
    std::cout << l2.isInvertible() << " inv\n";
    calibrated_jacobian_inverse=M.transpose().inverse();
    //std::cout << calibrated_jacobian_inverse << "\naaa" ;
    std::cout << "End of calib \n";
    
}

void DualOISCurveCalibration::PricingAfterCalibration (int n) {
    std::shared_ptr<aadc::AADCWorkSpace<mmType>> ws(aad_funcs.createWorkSpace());
    
    for (int i=0; i<calibrated_rates.size(); i++) {
        ws->val(calibrated_x_args[i])=aadc::mmSetConst<mmType>((calibrated_rates)[i]);
    }
    ScalarSetWSInputsFromVectorVisitor<mmType> visitor(*ws);
    portfolio_arg_prices.visitM(visitor, portfolio_init_prices);
    aad_funcs.forward(*ws);
    ws->resetDiff(); 
    ws->diff(portfolio_arg_res_prices.oisSwaps[n])=aadc::mmSetConst<mmType>(1.);
    aad_funcs.reverse(*ws);        
    VectorXd grad;
    grad.resize(calibrated_rates.size());

    for (int i=0; i<calibrated_x_args.size(); i++) {
        grad[i]=ws->diffp(calibrated_x_args[i])[0];
    }
    std::cout << " grad " << grad.size(); 
    std::cout << "Q_INd" << indices_cutted_Q;
    std::cout << "P_INd" << indices_cutted_P;
    grad=grad(indices_cutted_Q);
    std::cout << "\n grad converted " << grad.size();
    //std::cout << calibrated_jacobian.size() <<"\n";
    //std::cout << calibrated_jacobian->invertible() << "\n";
    //std::cout << (*calibrated_jacobian).inverse() << "\n";
    std::cout << grad.size() << " New grad \n";
    //std::cout << grad << "\n";
    std::cout <<"\n----------\n";
    //    std::cout << calibrated_jacobian_inverse << "\n";
    grad=calibrated_jacobian_inverse * grad;
    std::cout << grad << "\n";
    std::cout << "aaa\n";
    //std::cout << calibrated_jacobian_inverse;
}




void CurveCalibration::calibrateCurve()
{
    // optimisation code
    std::cout << "Start optmisation: " << std::endl;
    std::shared_ptr<aadc::AADCWorkSpace<mmType>> ws(aad_funcs.createWorkSpace());
    VectorXd x_vals;
    std::vector<aadc::AADCArgument> x_args;
    
    //
    std::vector<aadc::AADCResult> y_args;
    VectorXd y_vals;
    InitResults(y_vals, y_args);
    //
    
    InitParameters(x_vals, x_args);
    
    
    auto bfgs_func = [&] (const VectorXd& x, VectorXd& grad)  {
		for (int i=0; i<x.size(); i++) {
			ws->val(x_args[i])=aadc::mmSetConst<mmType>(x[i]);
		}
        ScalarSetWSInputsFromVectorVisitor<mmType> visitor(*ws);
        portfolio_arg_prices.visitM(visitor, portfolio_init_prices);
        //AADDebugStart();
		aad_funcs.forward(*ws);
	    for (int i=0; i<y_vals.size(); i++) {
            std::cout << i << " " <<ws->valp(y_args[i])[0]<<"\n";
    	}
    
    
    	ws->resetDiff(); 
		ws->diff(res_arg)=aadc::mmSetConst<mmType>(1.);
        aad_funcs.reverse(*ws);        
		//AADDebugStop();
        double *_f  = (ws->valp(res_arg));

		for (int i=0; i<x.size(); i++) {
			grad[i]=ws->diffp(x_args[i])[0];
            //std::cout << "diff " << ws->diff(x_args[i])[0] << "\n";
		}
	    std::cout << "\n fx= " <<  _f[0] << "; x = " << x.transpose() << std::endl;
		return _f[0];
	};

	LBFGSParam<double> param;
    param.epsilon = 1e-14;
    param.max_iterations = 1000;
    LBFGSSolver<double> solver(param);   
    double fx;
	int niter = solver.minimize(bfgs_func, x_vals, fx);
   	
	std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x_vals.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;
    std::cout << "x size:  = " << x_vals.size() << std::endl;
    for (int i=0; i<x_vals.size(); i++) {
        std::cout << x_vals[i] <<", "; 
    }
};