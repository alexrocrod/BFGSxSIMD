/*
 Copyright (C) 2017 Quaternion Risk Management Ltd
 Copyright (C) 2017 Aareal Bank AG

 All rights reserved.

 This file is part of ORE, a free-software/open-source library
 for transparent pricing and risk analysis - http://opensourcerisk.org

 ORE is free software: you can redistribute it and/or modify it
 under the terms of the Modified BSD License.  You should have received a
 copy of the license along with this program.
 The license is also available online at <http://opensourcerisk.org>

 This program is distributed on the basis that it will form a useful
 contribution to risk analytics and model standardisation, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE. See the license for more details.
*/

#include "ored/utilities/to_string.hpp"
#include <orea/app/reportwriter.hpp>

//FIXME: including all is slow and bad
#include <orea/orea.hpp>
#include <ored/ored.hpp>
#include <ored/portfolio/structuredtradeerror.hpp>
#include <ostream>
#include <ql/cashflows/averagebmacoupon.hpp>
#include <ql/cashflows/indexedcashflow.hpp>
#include <ql/cashflows/inflationcoupon.hpp>
#include <ql/errors.hpp>
#include <qle/cashflows/fxlinkedcashflow.hpp>
#include <stdio.h>

#include <cpprest/json.h>

using std::string;
using std::vector;
using QuantLib::Date;
using ore::data::to_string;

using Eigen::VectorXd;


namespace ore {
namespace analytics {

void ReportWriter::writeNpv(ore::data::Report& report, const std::string& baseCurrency,
                            boost::shared_ptr<Market> market, const std::string& configuration,
                            boost::shared_ptr<Portfolio> portfolio) {
    LOG("portfolio valuation");
    DayCounter dc = ActualActual();
    Date today = Settings::instance().evaluationDate();
    report.addColumn("TradeId", string())
        .addColumn("TradeType", string())
        .addColumn("Maturity", Date())
        .addColumn("MaturityTime", Real(), 6)
        .addColumn("NPV", Real(), 6)
        .addColumn("NpvCurrency", string())
        .addColumn("NPV(Base)", Real(), 6)
        .addColumn("BaseCurrency", string())
        .addColumn("Notional", Real(), 2)
        .addColumn("Notional(Base)", Real(), 2)
        .addColumn("NettingSet", string())
        .addColumn("CounterParty", string());



#ifdef AADC_QL
#define  AADC_QL_RUN
#endif

#ifdef AADC_QL_RUN
  
    std::vector<Real> outputs;  
    map<tuple<string, YieldCurveType, string>, AADCCurveData> curves_data;

#ifdef AADC_512
    //    typedef __m256d mmType;
        typedef __m512d mmType; 
#else
        typedef __m256d mmType;
#endif


    std::shared_ptr<aadc::AADCFunctions<mmType> > aad_funcsToday;

    if (true) 
    {
        aad_funcsToday = std::shared_ptr<aadc::AADCFunctions<mmType> >(
            new aadc::AADCFunctions<mmType>(
                {
                    {aadc::AADC_UseCompressedCode, 0}
                    , {aadc::AADC_CodeBufferSize, 256*1024}
                    , {aadc::AADC_NumCompressorThreads, 9}
                    //, {aadc::AADC_ReuseConstants, 0}
                    //, {aadc::AADC_InitInputDiff, 0}   // not reset input_diff
                }            
            )
        );
    }
    else   
        aad_funcsToday = std::make_shared<aadc::AADCFunctions<mmType>>(); // slower 


    tuple<string, YieldCurveType, string> key = std::make_tuple("xois_eur", YieldCurveType::Discount, "EUR"); // curve used

    std::vector<aadc::AADCResult> output_args;
    market->saveYieldCurvesData(curves_data,key);
    aad_funcsToday->startRecording();
    market->markYieldCurvesAsInput(curves_data,key);

    for (auto trade : portfolio->trades()) {
        string npvCcy = trade->npvCurrency();
        Real fx = 1.0;
        if (npvCcy != baseCurrency)
            fx = market->fxSpot(npvCcy + baseCurrency, configuration)->value();
        
        Real npv = trade->instrument()->NPV();
        Real out = npv*fx;
        output_args.push_back(out.markAsOutput());
    }
    aad_funcsToday->stopRecording();
    aad_funcsToday->printPassiveExtractLocations(std::cout, "ORE Today");
    aad_funcsToday->outStats(std::cout, "ORE Today");

    std::shared_ptr<aadc::AADCWorkSpace<mmType> > ws(aad_funcsToday->createWorkSpace());

    // Ci calculation with normal parametres
    std::cout << "\n\n Fill ws with values ater bootstrapping";

    AADCCurveData item = curves_data[key];

    int count = item.curve_args.size(); 
    std::cout << "\n\n" << std::get<0>(key) << " " << (int)std::get<1>(key) << " " << std::get<2>(key)<< "\n: ";
    
    for (size_t j = 0; j < item.curve_args.size(); j++) 
    {
        double num = item.curve_vals[j].val;
        std::cout << num << " ";
        ws->val(item.curve_args[j]) = aadc::mmSetConst<mmType>(num);
    }

    aad_funcsToday->forward(*ws);                     

    std::cout << "\n\n Initial Outputs: \n\n";

    double total_npv = 0.0;
    vector<double> aadc_npv(output_args.size(), 0.);
    std::cout.precision(10);

    for(size_t i = 0; i< output_args.size(); i++) 
    {
        aadc_npv[i] = ws->valp(output_args[i])[0];
        total_npv += aadc_npv[i];
        std::cout << aadc_npv[i] << "\n";
    }
    std::cout << "total_npv = " << total_npv <<"\n";
    auto out_json = web::json::value::object();

    std::vector<std::tuple<int,int,int>>  non_zero_deltas;  //gamma

    for(size_t i = 0; i < output_args.size(); i++) 
    {
        auto out_json_trade = web::json::value::object();
        out_json_trade["npv"] = aadc_npv[i];

        ws->resetDiff();
        ws->diff(output_args[i]) = aadc::mmSetConst<mmType>(1);

        aad_funcsToday->reverse(*ws); 

        
    #if 0  // no need to print all the data      

        vector<web::json::value> dict, dict_for_total;

        auto diff_json = web::json::value::object();

        for (auto& item : curves_data) 
        {
            double total_diff(0.0);
            vector<web::json::value> str;

            // std::cout << "\n\n" << std::get<0>(item.first) << " " << std::get<2>(item.first)<< "\n: ";

            for (size_t j = 0; j < item.second.curve_args.size(); j++) 
            {
                double val = ws->diff(item.second.curve_args[j])[0];

                if (val == 0) continue;

                //non_zero_deltas.push_back(std::make_tuple(i,ic,j));    //gamma

                total_diff += val;

                auto bucket = web::json::value::object(); 
                bucket["curve.type1"] = web::json::value(std::get<0>(item.first));
                bucket["curve.type2"] = web::json::value(int(std::get<1>(item.first)));
                bucket["curve.name"] = web::json::value(std::get<2>(item.first));
                bucket["time"] = web::json::value(item.second.curve_times[j].val);
                bucket["delta"]= web::json::value(val);
                dict.push_back(bucket);
            }
            if (total_diff == 0) continue;

            auto tot_bucket = web::json::value::object();

            //tot_bucket["curve.type"] = web::json::value(double(curves_iter_to_key[ic].first));
            //tot_bucket["curve.name"] = web::json::value(curves_iter_to_key[ic].second);
            tot_bucket["delta"]= web::json::value(total_diff);

            dict_for_total.push_back(tot_bucket);
        
        out_json_trade["sensitiv"]=web::json::value::array(dict);
        out_json_trade["total sensitiv"]=web::json::value::array(dict_for_total);
        }
    #endif
        out_json["TRADE_ID" + to_string(i)] = out_json_trade;
    }
    // std::cout << out_json;
    std::cout << "\n-----------------------------start test-----------------------------\n";

#ifdef AADC_512
    std::cout<<"using mm512d\n";
#else
    std::cout<<"using mm256d\n";
#endif
    // save initial parameters
    VectorXd params_init(count);

    for (size_t j = 0; j < item.curve_args.size(); j++)
        params_init[j] = item.curve_vals[j].val;
    
    VectorXd grad_rand = VectorXd::Random(count);

    // random x in [1 - delta_x, 1 + delta_x]
    double  delta_x = 0.1; // default: 0.1

    VectorXd x_bfgs = params_init * (VectorXd::Random(count) *  delta_x + VectorXd::Constant(count,1)); 

    const VectorXd x_perturb = x_bfgs;
    
    // BfgsOreOld<mmType> func(aad_funcsToday, ws, output_args, aadc_npv,  curves_data, key); // for traditional approach to LBFGS++
    BfgsOre<mmType> func(aad_funcsToday, ws, output_args, aadc_npv,  curves_data, key); // without debug information

    double fx_init = func(x_bfgs,grad_rand);

    x_bfgs = x_perturb; // reset x
    func.resetNcalls();

    double fx = 0;
    long time = 0;
    int niter = -1 , niter_LS = -1;
    int ncalls_R = -1, ncalls_F = -1;

    LBFGSpp::LBFGSParam<double> param;
    param.max_linesearch = 100;
    param.epsilon_rel = 3e-4; // 2e-4 for most, some cases need this to increase to 3e-4 or 5e-4

    LBFGSpp::LBFGSSolver<double, LBFGSpp::LineSearchBacktracking> solver(param);

    // func.setDelta(1e-9); // needed for BfgsOremmType, bracketing, mm256d and Wolfe or Strong Wolfe
    // func.setPolyFitOrder(0); // disable polyfit
    // func.setDgOrder(0); // disable finite differences
    
    auto start = std::chrono::high_resolution_clock::now();

    std::tie(niter, niter_LS) = solver.minimize(func, x_bfgs, fx); 

    time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count();
    
    ncalls_R = func.getNcallsR();
    ncalls_F = func.getNcallsF();

    std::cout.precision(17);
    std::cout << std::endl << "Initial function value: " << fx_init << std::endl << std::endl;
    std::cout.precision(6);
    std::cout << "LBFGS iterations: " << niter  << ", Line Search Iterations: " << niter_LS  << std::endl;
    std::cout << "Reverse Calls: "<<  ncalls_R  << ", Forward Calls: " << ncalls_F << std::endl;
    std::cout << "Time: " << time  << " microseconds, f(x) = " << fx << std::endl;
    std::cout << "Final point x = \n " << x_bfgs.transpose() << std::endl << std::endl;

#endif

    for (auto trade : portfolio->trades()) {
        string npvCcy = trade->npvCurrency();
        Real fx = 1.0;
        if (npvCcy != baseCurrency)
            fx = market->fxSpot(npvCcy + baseCurrency, configuration)->value();
        try {
            Real npv = trade->instrument()->NPV();
            QL_REQUIRE(std::isfinite(npv), "npv is not finite (" << npv << ")");
            Date maturity = trade->maturity();
            report.next()
                .add(trade->id())
                .add(trade->tradeType())
                .add(maturity)
                .add(maturity == QuantLib::Null<Date>() ? Null<Real>() : dc.yearFraction(today, maturity))
                .add(npv)
                .add(npvCcy)
                .add(npv * fx)
                .add(baseCurrency)
                .add(trade->notional())
                .add(trade->notional() * fx)
                .add(trade->envelope().nettingSetId())
                .add(trade->envelope().counterparty());
        } catch (std::exception& e) {
            ALOG(StructuredTradeErrorMessage(trade->id(), trade->tradeType(), "Error during trade pricing", e.what()));
            Date maturity = trade->maturity();
            report.next()
                .add(trade->id())
                .add(trade->tradeType())
                .add(maturity)
                .add(maturity == QuantLib::Null<Date>() ? Null<Real>() : dc.yearFraction(today, maturity))
                .add(Null<Real>())
                .add("#NA")
                .add(Null<Real>())
                .add("#NA")
                .add(Null<Real>())
                .add(Null<Real>())
                .add("#NA")
                .add("#NA");
        }
    }
    report.end();
    LOG("NPV file written");
}

void ReportWriter::writeCashflow(ore::data::Report& report, boost::shared_ptr<ore::data::Portfolio> portfolio,
                                 boost::shared_ptr<ore::data::Market> market, const std::string& configuration) {
    Date asof = Settings::instance().evaluationDate();
    bool write_discount_factor = market ? true : false;
    LOG("Writing cashflow report for " << asof);
    report.addColumn("TradeId", string())
        .addColumn("Type", string())
        .addColumn("CashflowNo", Size())
        .addColumn("LegNo", Size())
        .addColumn("PayDate", Date())
        .addColumn("FlowType", string())
        .addColumn("Amount", Real(), 4)
        .addColumn("Currency", string())
        .addColumn("Coupon", Real(), 10)
        .addColumn("Accrual", Real(), 10)
        .addColumn("fixingDate", Date())
        .addColumn("fixingValue", Real(), 10)
        .addColumn("Notional", Real(), 4);

    if (write_discount_factor) {
        report.addColumn("DiscountFactor", Real(), 10);
        report.addColumn("PresentValue", Real(), 10);
    }
    const vector<boost::shared_ptr<Trade>>& trades = portfolio->trades();

    for (Size k = 0; k < trades.size(); k++) {
        if (trades[k]->tradeType() == "Swaption" || trades[k]->tradeType() == "CapFloor") {
            WLOG("cashflow for " << trades[k]->tradeType() << " " << trades[k]->id() << " skipped");
            continue;
        }
        try {
            const vector<Leg>& legs = trades[k]->legs();
            for (size_t i = 0; i < legs.size(); i++) {
                const QuantLib::Leg& leg = legs[i];
                bool payer = trades[k]->legPayers()[i];
                string ccy = trades[k]->legCurrencies()[i];
                Handle<YieldTermStructure> discountCurve;
                if (write_discount_factor)
                    discountCurve = market->discountCurve(ccy, configuration);
                for (size_t j = 0; j < leg.size(); j++) {
                    boost::shared_ptr<QuantLib::CashFlow> ptrFlow = leg[j];
                    Date payDate = ptrFlow->date();
                    if (payDate >= asof) {
                        Real amount = ptrFlow->amount();
                        string flowType = "";
                        if (payer)
                            amount *= -1.0;
                        std::string ccy = trades[k]->legCurrencies()[i];
                        boost::shared_ptr<QuantLib::Coupon> ptrCoupon =
                            boost::dynamic_pointer_cast<QuantLib::Coupon>(ptrFlow);
                        Real coupon;
                        Real accrual;
                        Real notional;
                        if (ptrCoupon) {
                            coupon = ptrCoupon->rate();
                            accrual = ptrCoupon->accrualPeriod();
                            notional = ptrCoupon->nominal();
                            flowType = "Interest";
                        } else {
                            coupon = Null<Real>();
                            accrual = Null<Real>();
                            notional = Null<Real>();
                            flowType = "Notional";
                        }
                        // This BMA part here (and below) is necessary because the fixingDay() method of
                        // AverageBMACoupon returns an exception rather than the last fixing day of the period.
                        boost::shared_ptr<AverageBMACoupon> ptrBMA =
                            boost::dynamic_pointer_cast<QuantLib::AverageBMACoupon>(ptrFlow);
                        boost::shared_ptr<QuantLib::FloatingRateCoupon> ptrFloat =
                            boost::dynamic_pointer_cast<QuantLib::FloatingRateCoupon>(ptrFlow);
                        boost::shared_ptr<QuantLib::InflationCoupon> ptrInfl =
                            boost::dynamic_pointer_cast<QuantLib::InflationCoupon>(ptrFlow);
                        boost::shared_ptr<QuantLib::IndexedCashFlow> ptrIndCf =
                            boost::dynamic_pointer_cast<QuantLib::IndexedCashFlow>(ptrFlow);
                        boost::shared_ptr<QuantExt::FXLinkedCashFlow> ptrFxlCf =
                            boost::dynamic_pointer_cast<QuantExt::FXLinkedCashFlow>(ptrFlow);
                        Date fixingDate;
                        Real fixingValue;
                        if (ptrBMA) {
                            // We return the last fixing inside the coupon period
                            fixingDate = ptrBMA->fixingDates().end()[-2];
                            fixingValue = ptrBMA->pricer()->swapletRate();
                            if (ptrBMA->index()->pastFixing(fixingDate) == Null<Real>())
                                flowType = "BMAaverage";
                        } else if (ptrFloat) {
                            fixingDate = ptrFloat->fixingDate();
                            fixingValue = ptrFloat->index()->fixing(fixingDate);
                            if (ptrFloat->index()->pastFixing(fixingDate) == Null<Real>())
                                flowType = "InterestProjected";
                        } else if (ptrInfl) {
                            fixingDate = ptrInfl->fixingDate();
                            fixingValue = ptrInfl->index()->fixing(fixingDate);
                            flowType = "Inflation";
                        } else if (ptrIndCf) {
                            fixingDate = ptrIndCf->fixingDate();
                            fixingValue = ptrIndCf->index()->fixing(fixingDate);
                            flowType = "Index";
                        } else if (ptrFxlCf) {
                            fixingDate = ptrFxlCf->fxFixingDate();
                            fixingValue = ptrFxlCf->fxRate();
                        } else {
                            fixingDate = Null<Date>();
                            fixingValue = Null<Real>();
                        }
                        report.next()
                            .add(trades[k]->id())
                            .add(trades[k]->tradeType())
                            .add(j + 1)
                            .add(i)
                            .add(payDate)
                            .add(flowType)
                            .add(amount)
                            .add(ccy)
                            .add(coupon)
                            .add(accrual)
                            .add(fixingDate)
                            .add(fixingValue)
                            .add(notional);

                        if (write_discount_factor) {
                            Real discountFactor = discountCurve->discount(payDate);
                            report.add(discountFactor);
                            Real presentValue = discountFactor * amount;
                            report.add(presentValue);
                        }
                    }
                }
            }
        } catch (std::exception& e) {
            LOG("Exception writing cashflow report : " << e.what());
        } catch (...) {
            LOG("Exception writing cashflow report : Unkown Exception");
        }
    }
    report.end();
    LOG("Cashflow report written");
}

void ReportWriter::writeCurves(ore::data::Report& report, const std::string& configID, const DateGrid& grid,
                               const TodaysMarketParameters& marketConfig, const boost::shared_ptr<Market>& market,
                               const bool continueOnError) {
    LOG("Write curves... ");

    QL_REQUIRE(marketConfig.hasConfiguration(configID), "curve configuration " << configID << " not found");

    map<string, string> discountCurves = marketConfig.mapping(MarketObject::DiscountCurve, configID);
    map<string, string> YieldCurves = marketConfig.mapping(MarketObject::YieldCurve, configID);
    map<string, string> indexCurves = marketConfig.mapping(MarketObject::IndexCurve, configID);
    map<string, string> zeroInflationIndices, defaultCurves;
    if (marketConfig.hasMarketObject(MarketObject::ZeroInflationCurve))
        zeroInflationIndices = marketConfig.mapping(MarketObject::ZeroInflationCurve, configID);
    if (marketConfig.hasMarketObject(MarketObject::DefaultCurve))
        defaultCurves = marketConfig.mapping(MarketObject::DefaultCurve, configID);

    vector<Handle<YieldTermStructure>> yieldCurves;
    vector<Handle<ZeroInflationIndex>> zeroInflationFixings;
    vector<Handle<DefaultProbabilityTermStructure>> probabilityCurves;

    report.addColumn("Tenor", Period()).addColumn("Date", Date());

    for (auto it : discountCurves) {
        DLOG("discount curve - " << it.first);
        try {
            yieldCurves.push_back(market->discountCurve(it.first, configID));
            report.addColumn(it.first, Real(), 15);
        } catch (const std::exception& e) {
            if (continueOnError) {
                WLOG("skip this curve: " << e.what());
            } else {
                QL_FAIL(e.what());
            }
        }
    }
    for (auto it : YieldCurves) {
        DLOG("yield curve - " << it.first);
        try {
            yieldCurves.push_back(market->yieldCurve(it.first, configID));
            report.addColumn(it.first, Real(), 15);
        } catch (const std::exception& e) {
            if (continueOnError) {
                WLOG("skip this curve: " << e.what());
            } else {
                QL_FAIL(e.what());
            }
        }
    }
    for (auto it : indexCurves) {
        DLOG("index curve - " << it.first);
        try {
            yieldCurves.push_back(market->iborIndex(it.first, configID)->forwardingTermStructure());
            report.addColumn(it.first, Real(), 15);
        } catch (const std::exception& e) {
            if (continueOnError) {
                WLOG("skip this curve: " << e.what());
            } else {
                QL_FAIL(e.what());
            }
        }
    }
    for (auto it : zeroInflationIndices) {
        DLOG("inflation curve - " << it.first);
        try {
            zeroInflationFixings.push_back(market->zeroInflationIndex(it.first, configID));
            report.addColumn(it.first, Real(), 15);
        } catch (const std::exception& e) {
            if (continueOnError) {
                WLOG("skip this curve: " << e.what());
            } else {
                QL_FAIL(e.what());
            }
        }
    }
    for (auto it : defaultCurves) {
        DLOG("default curve - " << it.first);
        try {
            probabilityCurves.push_back(market->defaultCurve(it.first, configID));
            report.addColumn(it.first, Real(), 15);
        } catch (const std::exception& e) {
            if (continueOnError) {
                WLOG("skip this curve: " << e.what());
            } else {
                QL_FAIL(e.what());
            }
        }
    }

    for (Size j = 0; j < grid.size(); ++j) {
        Date date = grid[j];
        report.next().add(grid.tenors()[j]).add(date);
        for (Size i = 0; i < yieldCurves.size(); ++i)
            report.add(yieldCurves[i]->discount(date));
        for (Size i = 0; i < zeroInflationFixings.size(); ++i)
            report.add(zeroInflationFixings[i]->fixing(date));
        for (Size i = 0; i < probabilityCurves.size(); ++i)
            report.add(probabilityCurves[i]->survivalProbability(date));
    }
    report.end();
}

void ReportWriter::writeTradeExposures(ore::data::Report& report, boost::shared_ptr<PostProcess> postProcess,
                                       const string& tradeId) {
    const vector<Date> dates = postProcess->cube()->dates();
    Date today = Settings::instance().evaluationDate();
    DayCounter dc = ActualActual();
    const vector<Real>& epe = postProcess->tradeEPE(tradeId);
    const vector<Real>& ene = postProcess->tradeENE(tradeId);
    const vector<Real>& ee_b = postProcess->tradeEE_B(tradeId);
    const vector<Real>& eee_b = postProcess->tradeEEE_B(tradeId);
    const vector<Real>& pfe = postProcess->tradePFE(tradeId);
    const vector<Real>& aepe = postProcess->allocatedTradeEPE(tradeId);
    const vector<Real>& aene = postProcess->allocatedTradeENE(tradeId);
    report.addColumn("TradeId", string())
        .addColumn("Date", Date())
        .addColumn("Time", Real(), 6)
        .addColumn("EPE", Real())
        .addColumn("ENE", Real())
        .addColumn("AllocatedEPE", Real())
        .addColumn("AllocatedENE", Real())
        .addColumn("PFE", Real())
        .addColumn("BaselEE", Real())
        .addColumn("BaselEEE", Real());

    report.next()
        .add(tradeId)
        .add(today)
        .add(Real(0.0))
        .add(epe[0])
        .add(ene[0])
        .add(aepe[0])
        .add(aene[0])
        .add(pfe[0])
        .add(ee_b[0])
        .add(eee_b[0]);
    for (Size j = 0; j < dates.size(); ++j) {
        Time time = dc.yearFraction(today, dates[j]);
        report.next()
            .add(tradeId)
            .add(dates[j])
            .add(time)
            .add(epe[j + 1])
            .add(ene[j + 1])
            .add(aepe[j + 1])
            .add(aene[j + 1])
            .add(pfe[j + 1])
            .add(ee_b[j + 1])
            .add(eee_b[j + 1]);
    }
    report.end();
}

void ReportWriter::writeNettingSetExposures(ore::data::Report& report, boost::shared_ptr<PostProcess> postProcess,
                                            const string& nettingSetId) {
    const vector<Date> dates = postProcess->cube()->dates();
    Date today = Settings::instance().evaluationDate();
    DayCounter dc = ActualActual();
    const vector<Real>& epe = postProcess->netEPE(nettingSetId);
    const vector<Real>& ene = postProcess->netENE(nettingSetId);
    const vector<Real>& ee_b = postProcess->netEE_B(nettingSetId);
    const vector<Real>& eee_b = postProcess->netEEE_B(nettingSetId);
    const vector<Real>& pfe = postProcess->netPFE(nettingSetId);
    const vector<Real>& ecb = postProcess->expectedCollateral(nettingSetId);
    report.addColumn("NettingSet", string())
        .addColumn("Date", Date())
        .addColumn("Time", Real(), 6)
        .addColumn("EPE", Real(), 2)
        .addColumn("ENE", Real(), 2)
        .addColumn("PFE", Real(), 2)
        .addColumn("ExpectedCollateral", Real(), 2)
        .addColumn("BaselEE", Real(), 2)
        .addColumn("BaselEEE", Real(), 2);

    report.next()
        .add(nettingSetId)
        .add(today)
        .add(Real(0.0))
        .add(epe[0])
        .add(ene[0])
        .add(pfe[0])
        .add(ecb[0])
        .add(ee_b[0])
        .add(eee_b[0]);
    for (Size j = 0; j < dates.size(); ++j) {
        Real time = dc.yearFraction(today, dates[j]);
        report.next()
            .add(nettingSetId)
            .add(dates[j])
            .add(time)
            .add(epe[j + 1])
            .add(ene[j + 1])
            .add(pfe[j + 1])
            .add(ecb[j + 1])
            .add(ee_b[j + 1])
            .add(eee_b[j + 1]);
    }
    report.end();
}

void ReportWriter::writeXVA(ore::data::Report& report, const string& allocationMethod,
                            boost::shared_ptr<Portfolio> portfolio, boost::shared_ptr<PostProcess> postProcess) {
    const vector<Date> dates = postProcess->cube()->dates();
    DayCounter dc = ActualActual();
    report.addColumn("TradeId", string())
        .addColumn("NettingSetId", string())
        .addColumn("CVA", Real(), 2)
        .addColumn("DVA", Real(), 2)
        .addColumn("FBA", Real(), 2)
        .addColumn("FCA", Real(), 2)
        .addColumn("FBAexOwnSP", Real(), 2)
        .addColumn("FCAexOwnSP", Real(), 2)
        .addColumn("FBAexAllSP", Real(), 2)
        .addColumn("FCAexAllSP", Real(), 2)
        .addColumn("COLVA", Real(), 2)
        .addColumn("MVA", Real(), 2)
        .addColumn("OurKVACCR", Real(), 2)
        .addColumn("TheirKVACCR", Real(), 2)
        .addColumn("OurKVACVA", Real(), 2)
        .addColumn("TheirKVACVA", Real(), 2)
        .addColumn("CollateralFloor", Real(), 2)
        .addColumn("AllocatedCVA", Real(), 2)
        .addColumn("AllocatedDVA", Real(), 2)
        .addColumn("AllocationMethod", string())
        .addColumn("BaselEPE", Real(), 2)
        .addColumn("BaselEEPE", Real(), 2);

    for (auto n : postProcess->nettingSetIds()) {
        report.next()
            .add("")
            .add(n)
            .add(postProcess->nettingSetCVA(n))
            .add(postProcess->nettingSetDVA(n))
            .add(postProcess->nettingSetFBA(n))
            .add(postProcess->nettingSetFCA(n))
            .add(postProcess->nettingSetFBA_exOwnSP(n))
            .add(postProcess->nettingSetFCA_exOwnSP(n))
            .add(postProcess->nettingSetFBA_exAllSP(n))
            .add(postProcess->nettingSetFCA_exAllSP(n))
            .add(postProcess->nettingSetCOLVA(n))
            .add(postProcess->nettingSetMVA(n))
            .add(postProcess->nettingSetOurKVACCR(n))
            .add(postProcess->nettingSetTheirKVACCR(n))
            .add(postProcess->nettingSetOurKVACVA(n))
            .add(postProcess->nettingSetTheirKVACVA(n))
            .add(postProcess->nettingSetCollateralFloor(n))
            .add(postProcess->nettingSetCVA(n))
            .add(postProcess->nettingSetDVA(n))
            .add(allocationMethod)
            .add(postProcess->netEPE_B(n))
            .add(postProcess->netEEPE_B(n));

        for (Size k = 0; k < portfolio->trades().size(); ++k) {
            string tid = portfolio->trades()[k]->id();
            string nid = portfolio->trades()[k]->envelope().nettingSetId();
            if (nid != n)
                continue;
            report.next()
                .add(tid)
                .add(nid)
                .add(postProcess->tradeCVA(tid))
                .add(postProcess->tradeDVA(tid))
                .add(postProcess->tradeFBA(tid))
                .add(postProcess->tradeFCA(tid))
                .add(postProcess->tradeFBA_exOwnSP(tid))
                .add(postProcess->tradeFCA_exOwnSP(tid))
                .add(postProcess->tradeFBA_exAllSP(tid))
                .add(postProcess->tradeFCA_exAllSP(tid))
                .add(Null<Real>())
                .add(Null<Real>())
                .add(Null<Real>())
                .add(Null<Real>())
                .add(Null<Real>())
                .add(Null<Real>())
                .add(Null<Real>())
                .add(postProcess->allocatedTradeCVA(tid))
                .add(postProcess->allocatedTradeDVA(tid))
                .add(allocationMethod)
                .add(postProcess->tradeEPE_B(tid))
                .add(postProcess->tradeEEPE_B(tid));
        }
    }
    report.end();
}

void ReportWriter::writeNettingSetColva(ore::data::Report& report, boost::shared_ptr<PostProcess> postProcess,
                                        const string& nettingSetId) {
    const vector<Date> dates = postProcess->cube()->dates();
    Date today = Settings::instance().evaluationDate();
    DayCounter dc = ActualActual();
    const vector<Real>& collateral = postProcess->expectedCollateral(nettingSetId);
    const vector<Real>& colvaInc = postProcess->colvaIncrements(nettingSetId);
    const vector<Real>& floorInc = postProcess->collateralFloorIncrements(nettingSetId);
    Real colva = postProcess->nettingSetCOLVA(nettingSetId);
    Real floorValue = postProcess->nettingSetCollateralFloor(nettingSetId);

    report.addColumn("NettingSet", string())
        .addColumn("Date", Date())
        .addColumn("Time", Real(), 4)
        .addColumn("CollateralBalance", Real(), 4)
        .addColumn("COLVA Increment", Real(), 4)
        .addColumn("COLVA", Real(), 4)
        .addColumn("CollateralFloor Increment", Real(), 4)
        .addColumn("CollateralFloor", Real(), 4);

    report.next()
        .add(nettingSetId)
        .add(Null<Date>())
        .add(Null<Real>())
        .add(Null<Real>())
        .add(Null<Real>())
        .add(colva)
        .add(Null<Real>())
        .add(floorValue);
    Real colvaSum = 0.0;
    Real floorSum = 0.0;
    for (Size j = 0; j < dates.size(); ++j) {
        Real time = dc.yearFraction(today, dates[j]);
        colvaSum += colvaInc[j + 1];
        floorSum += floorInc[j + 1];
        report.next()
            .add(nettingSetId)
            .add(dates[j])
            .add(time)
            .add(collateral[j + 1])
            .add(colvaInc[j + 1])
            .add(colvaSum)
            .add(floorInc[j + 1])
            .add(floorSum);
    }
    report.end();
}

void ReportWriter::writeAggregationScenarioData(ore::data::Report& report, const AggregationScenarioData& data) {
    report.addColumn("Date", Size()).addColumn("Scenario", Size());
    for (auto const& k : data.keys()) {
        std::string tmp = ore::data::to_string(k.first) + k.second;
        report.addColumn(tmp.c_str(), Real(), 8);
    }
    for (Size d = 0; d < data.dimDates(); ++d) {
        for (Size s = 0; s < data.dimSamples(); ++s) {
            report.next();
            report.add(d).add(s);
            for (auto const& k : data.keys()) {
                report.add(data.get(d, s, k.first, k.second));
            }
        }
    }
    report.end();
}

void ReportWriter::writeScenarioReport(ore::data::Report& report, const boost::shared_ptr<SensitivityCube>& sensitivityCube,
                                       Real outputThreshold) {

    LOG("Writing Scenario report");

    report.addColumn("TradeId", string());
    report.addColumn("Factor", string());
    report.addColumn("Up/Down", string());
    report.addColumn("Base NPV", Real(), 2);
    report.addColumn("Scenario NPV", Real(), 2);
    report.addColumn("Difference", Real(), 2);

    auto scenarioDescriptions = sensitivityCube->scenarioDescriptions();
    auto tradeIds = sensitivityCube->tradeIds();
    auto npvCube = sensitivityCube->npvCube();

    for (Size i = 0; i < tradeIds.size(); i++) {
        Real baseNpv = npvCube->getT0(i);
        auto tradeId = tradeIds[i];

        for (Size j = 0; j < scenarioDescriptions.size(); j++) {
            auto scenarioDescription = scenarioDescriptions[j];

            Real scenarioNpv = npvCube->get(i, j);
            Real difference = scenarioNpv - baseNpv;

            if (std::fabs(difference) > outputThreshold) {
                report.next();
                report.add(tradeId);
                report.add(scenarioDescription.factors());
                report.add(scenarioDescription.typeString());
                report.add(baseNpv);
                report.add(scenarioNpv);
                report.add(difference);
            }
            else if (!std::isfinite(difference)) {
                // TODO: is this needed?
                ALOG("sensitivity scenario for trade " << tradeId << ", factor " << scenarioDescription.factors()
                    << " is not finite (" << difference << ")");
            }
        }
    }

    report.end();
    LOG("Scenario report finished");
}

void ReportWriter::writeSensitivityReport(Report& report, const boost::shared_ptr<SensitivityStream>& ss,
                                          Real outputThreshold) {

    LOG("Writing Sensitivity report");

    report.addColumn("TradeId", string());
    report.addColumn("IsPar", string());
    report.addColumn("Factor_1", string());
    report.addColumn("ShiftSize_1", Real(), 6);
    report.addColumn("Factor_2", string());
    report.addColumn("ShiftSize_2", Real(), 6);
    report.addColumn("Currency", string());
    report.addColumn("Base NPV", Real(), 2);
    report.addColumn("Delta", Real(), 2);
    report.addColumn("Gamma", Real(), 2);

    // Make sure that we are starting from the start
    ss->reset();
    while (SensitivityRecord sr = ss->next()) {
        if (std::fabs(sr.delta) > outputThreshold || std::fabs(sr.gamma) > outputThreshold) {
            report.next();
            report.add(sr.tradeId);
            report.add(to_string(sr.isPar));
            report.add(reconstructFactor(sr.key_1, sr.desc_1));
            report.add(sr.shift_1);
            report.add(reconstructFactor(sr.key_2, sr.desc_2));
            report.add(sr.shift_2);
            report.add(sr.currency);
            report.add(sr.baseNpv);
            report.add(sr.delta);
            report.add(sr.gamma);
        } else if (!std::isfinite(sr.delta) || !std::isfinite(sr.gamma)) {
            // TODO: Again, is this needed?
            ALOG("sensitivity record has infinite values: " << sr);
        }
    }

    report.end();
    LOG("Sensitivity report finished");
}

} // namespace analytics
} // namespace ore
