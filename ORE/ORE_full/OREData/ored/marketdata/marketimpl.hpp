/*
 Copyright (C) 2016 Quaternion Risk Management Ltd
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

/*! \file ored/marketdata/marketimpl.hpp
    \brief An implementation of the Market class that stores the required objects in maps
    \ingroup marketdata
*/

#pragma once

#include "ored/utilities/xmlutils.hpp"
#include <ored/configuration/conventions.hpp>
#include <ored/marketdata/fxtriangulation.hpp>
#include <ored/marketdata/market.hpp>

#include <qle/indexes/inflationindexobserver.hpp>

#include <map>

#ifdef AADC_QL
#include <aadc/aadc.h>
#include <aadc/aadc_thread_sync.h>
#include <aadc/aadc_debug.h>

#endif

namespace ore {
namespace data {
using namespace QuantLib;
using ore::data::Convention;
using ore::data::Conventions;
using std::string;
using std::map;
using std::pair;
using std::tuple;

// TODO: rename class
//! Market Implementation
/*!
  The MarketImpl class differs from the Market base class in that it contains concrete maps
  of term structures, and it implements the interface.

  \ingroup marketdata
 */
class MarketImpl : public Market {
public:
    //! Default constructor
    MarketImpl() {}
    MarketImpl(const Conventions& conventions) : conventions_(conventions) {
        // if no fx spots are defined we still need an empty triangulation
        fxSpots_[Market::defaultConfiguration] = FXTriangulation();
    }

    //! \name Market interface
    //@{
    //! Get the asof Date
    Date asofDate() const { return asof_; }

    //! Yield Curves
    Handle<YieldTermStructure> yieldCurve(const YieldCurveType& type, const string& ccy,
                                          const string& configuration = Market::defaultConfiguration) const;
    Handle<YieldTermStructure> discountCurve(const string& ccy,
                                             const string& configuration = Market::defaultConfiguration) const;
    Handle<YieldTermStructure> yieldCurve(const string& name,
                                          const string& configuration = Market::defaultConfiguration) const;
    Handle<IborIndex> iborIndex(const string& indexName,
                                const string& configuration = Market::defaultConfiguration) const;
    Handle<SwapIndex> swapIndex(const string& indexName,
                                const string& configuration = Market::defaultConfiguration) const;

    //! Swaptions
    Handle<QuantLib::SwaptionVolatilityStructure>
    swaptionVol(const string& ccy, const string& configuration = Market::defaultConfiguration) const;
    const string shortSwapIndexBase(const string& ccy,
                                    const string& configuration = Market::defaultConfiguration) const;
    const string swapIndexBase(const string& ccy, const string& configuration = Market::defaultConfiguration) const;

    //! Yield volatility
    Handle<QuantLib::SwaptionVolatilityStructure>
        yieldVol(const string& securityID, const string& configuration = Market::defaultConfiguration) const;

    //! FX
    Handle<Quote> fxSpot(const string& ccypair, const string& configuration = Market::defaultConfiguration) const;
    Handle<BlackVolTermStructure> fxVol(const string& ccypair,
                                        const string& configuration = Market::defaultConfiguration) const;

    //! Default Curves and Recovery Rates
    Handle<DefaultProbabilityTermStructure>
    defaultCurve(const string&, const string& configuration = Market::defaultConfiguration) const;
    Handle<Quote> recoveryRate(const string&, const string& configuration = Market::defaultConfiguration) const;

    //! CDS volatilities
    Handle<BlackVolTermStructure> cdsVol(const string& name,
                                         const string& configuration = Market::defaultConfiguration) const;

    //! Base correlation structures
    Handle<BaseCorrelationTermStructure<BilinearInterpolation>>
    baseCorrelation(const string& name, const string& configuration = Market::defaultConfiguration) const;

    //! CapFloor volatilities
    Handle<OptionletVolatilityStructure> capFloorVol(const string& ccy,
                                                     const string& configuration = Market::defaultConfiguration) const;

    //! YoY Inflation CapFloor volatilities
    Handle<QuantExt::YoYOptionletVolatilitySurface>
    yoyCapFloorVol(const string& name, const string& configuration = Market::defaultConfiguration) const;

    //! Inflation Indexes
    virtual Handle<ZeroInflationIndex>
    zeroInflationIndex(const string& indexName, const string& configuration = Market::defaultConfiguration) const;
    virtual Handle<YoYInflationIndex>
    yoyInflationIndex(const string& indexName, const string& configuration = Market::defaultConfiguration) const;

    //! CPI Inflation Cap Floor Price Surfaces
    virtual Handle<CPICapFloorTermPriceSurface>
    cpiInflationCapFloorPriceSurface(const string& indexName,
                                     const string& configuration = Market::defaultConfiguration) const;

    //! Inflation Cap Floor Volatility Surfaces
    virtual Handle<CPIVolatilitySurface>
    cpiInflationCapFloorVolatilitySurface(const string& indexName,
                                          const string& configuration = Market::defaultConfiguration) const;

    //! YoY Inflation Cap Floor Price Surfaces
    virtual Handle<YoYCapFloorTermPriceSurface>
    yoyInflationCapFloorPriceSurface(const string& indexName,
                                     const string& configuration = Market::defaultConfiguration) const;

    //! Equity curves
    Handle<Quote> equitySpot(const string& eqName, const string& configuration = Market::defaultConfiguration) const;
    Handle<QuantExt::EquityIndex> equityCurve(const string& eqName,
                                              const string& configuration = Market::defaultConfiguration) const;

    Handle<YieldTermStructure> equityDividendCurve(const string& eqName,
                                                   const string& configuration = Market::defaultConfiguration) const;

    //! Equity volatilities
    Handle<BlackVolTermStructure> equityVol(const string& eqName,
                                            const string& configuration = Market::defaultConfiguration) const;

    //! Equity forecasting curves
    Handle<YieldTermStructure> equityForecastCurve(const string& eqName,
                                                   const string& configuration = Market::defaultConfiguration) const;

    //! Bond Spreads
    Handle<Quote> securitySpread(const string& securityID,
                                 const string& configuration = Market::defaultConfiguration) const;

    //! Cpi Base Quotes
    Handle<QuantExt::InflationIndexObserver> baseCpis(const string& index,
                                                      const string& configuration = Market::defaultConfiguration) const;

    //! Commodity curves
    QuantLib::Handle<QuantLib::Quote> commoditySpot(const string& commodityName,
                                                    const string& configuration = Market::defaultConfiguration) const;

    QuantLib::Handle<QuantExt::PriceTermStructure>
    commodityPriceCurve(const string& commodityName, const string& configuration = Market::defaultConfiguration) const;

    QuantLib::Handle<QuantLib::BlackVolTermStructure>
    commodityVolatility(const string& commodityName, const string& configuration = Market::defaultConfiguration) const;
    //@}

    //! Correlation curves
    Handle<QuantExt::CorrelationTermStructure>
    correlationCurve(const string& index1, const string& index2,
                     const string& configuration = Market::defaultConfiguration) const;
    //! \name Conditional Prepayment Rates
    //@{
    QuantLib::Handle<Quote> cpr(const string& securityID,
                                const string& configuration = Market::defaultConfiguration) const;
    //@}

    //! \name Disable copying
    //@{
    MarketImpl(const MarketImpl&) = delete;
    MarketImpl& operator=(const MarketImpl&) = delete;
    //@}

    //! Send an explicit update() call to all term structures
    void refresh(const string& configuration = Market::defaultConfiguration);

    // void markYieldCurvesAsInput (map<tuple<string, YieldCurveType, string>, AADCCurveData>& curves_data,bool full_curves=true) {  
    void markYieldCurvesAsInput (map<tuple<string, YieldCurveType, string>, AADCCurveData>& curves_data) {  
        auto key = std::make_tuple("xois_eur", YieldCurveType::Discount, "EUR"); 
        for (auto& item: yieldCurves_) {
            // if (full_curves || item.first==key){
            if (item.first==key){
                // negation of the condition in the next string is (== Or ==), i.e. if (== Or ==), then mark this curve. To exclude equities.
                if (std::get<1>(item.first) != YieldCurveType::Discount && std::get<1>(item.first) != YieldCurveType::Yield) continue; 
                //|| std::get<0>(item.first) != "default") continue;
                aadc::VectorArg vec_args;
                item.second->markVecAsInput(vec_args);
                curves_data[item.first].curve_args=vec_args;
            }
        }
    }
    void markYieldCurvesAsInput (map<tuple<string, YieldCurveType, string>, AADCCurveData>& curves_data, tuple<string, YieldCurveType, string> &key) {  
        for (auto& item: yieldCurves_) {
            if (item.first==key){
                // negation of the condition in the next string is (== Or ==), i.e. if (== Or ==), then mark this curve. To exclude equities.
                if (std::get<1>(item.first) != YieldCurveType::Discount && std::get<1>(item.first) != YieldCurveType::Yield) continue; 
                aadc::VectorArg vec_args;
                item.second->markVecAsInput(vec_args);
                curves_data[item.first].curve_args=vec_args;
            }
        }
    }

    void markYieldCurvesOutput(map<tuple<string, YieldCurveType, string>, aadc::VectorRes>& curves_outargs_) {
        for (auto& item: yieldCurves_) {
            // negation of the condition in the next string is (== Or ==), i.e. if (== Or ==), then mark this curve. To exclude equities.
            if (std::get<1>(item.first) != YieldCurveType::Discount && std::get<1>(item.first) != YieldCurveType::Yield) continue; 
            //|| std::get<0>(item.first) != "default") continue;
            aadc::VectorRes vec_args;
            item.second->markVecAsOutput(vec_args);
            curves_outargs_[item.first] = vec_args;
        }
    }

    //void saveYieldCurvesData (map<tuple<string, YieldCurveType, string>, AADCCurveData>& curves_data,bool full_curves=true) {
    void saveYieldCurvesData (map<tuple<string, YieldCurveType, string>, AADCCurveData>& curves_data) {
        auto key = std::make_tuple("xois_eur", YieldCurveType::Discount, "EUR"); 
        for (auto& item: yieldCurves_) {
            // if (full_curves || item.first==key){
            if (item.first==key){
                if (std::get<1>(item.first) != YieldCurveType::Discount  && std::get<1>(item.first) != YieldCurveType::Yield ) continue;
                vector<Real> vec_vals(item.second->data());
                AADCCurveData& tmp= curves_data[item.first]; 
                // for (int i=1; i< vec_vals.size(); i++) {
                //     //vec_vals[i]=-std::log(vec_vals[i]) / item.second->times()[i]; // No transition to zero-rates 
                //                                                                     // in original curves of TodaysMarket
                // }
                tmp.curve_vals=vec_vals;
                tmp.curve_times=item.second->times();
            }
        }
    }

    void saveYieldCurvesData (map<tuple<string, YieldCurveType, string>, AADCCurveData>& curves_data, tuple<string, YieldCurveType, string> &key) {
        for (auto& item: yieldCurves_) {
            if (item.first==key){
                if (std::get<1>(item.first) != YieldCurveType::Discount  && std::get<1>(item.first) != YieldCurveType::Yield ) continue;
                vector<Real> vec_vals(item.second->data());
                AADCCurveData& tmp= curves_data[item.first]; 
                tmp.curve_vals=vec_vals;
                tmp.curve_times=item.second->times();
            }
        }
    }


    void shift (int i, Real coeff, tuple<string, YieldCurveType, string> curve_key) {
        
        yieldCurves_[curve_key]->data()[yieldCurves_[curve_key]->data().size()-i-1]+=coeff;
        
    
        yieldCurves_[curve_key]->data()[yieldCurves_[curve_key]->data().size()-i-1]= 
            yieldCurves_[curve_key]->data()[yieldCurves_[curve_key]->data().size()-i-1];
       // yieldCurves_[curve_key]->data()[yieldCurves_[curve_key]->data().size()-i-1]=AADC_PRINT(
        //    yieldCurves_[curve_key]->data()[yieldCurves_[curve_key]->data().size()-i-1]
        //);
        //std::cout << MyDebugPrint::call_counter << "\n";
        yieldCurves_[curve_key]->CurveUpdate();
    }

protected:
    Date asof_;
    // maps (configuration, key) => term structure
    map<tuple<string, YieldCurveType, string>, Handle<YieldTermStructure>> yieldCurves_;
    map<pair<string, string>, Handle<IborIndex>> iborIndices_;
    map<pair<string, string>, Handle<SwapIndex>> swapIndices_;
    map<pair<string, string>, Handle<QuantLib::SwaptionVolatilityStructure>> swaptionCurves_;
    map<pair<string, string>, pair<string, string>> swaptionIndexBases_;
    map<pair<string, string>, Handle<QuantLib::SwaptionVolatilityStructure>> yieldVolCurves_;
    map<string, FXTriangulation> fxSpots_;
    mutable map<pair<string, string>, Handle<BlackVolTermStructure>> fxVols_;
    map<pair<string, string>, Handle<DefaultProbabilityTermStructure>> defaultCurves_;
    map<pair<string, string>, Handle<BlackVolTermStructure>> cdsVols_;
    map<pair<string, string>, Handle<BaseCorrelationTermStructure<BilinearInterpolation>>> baseCorrelations_;
    map<pair<string, string>, Handle<Quote>> recoveryRates_;
    map<pair<string, string>, Handle<OptionletVolatilityStructure>> capFloorCurves_;
    map<pair<string, string>, Handle<QuantExt::YoYOptionletVolatilitySurface>> yoyCapFloorVolSurfaces_;
    map<pair<string, string>, Handle<ZeroInflationIndex>> zeroInflationIndices_;
    map<pair<string, string>, Handle<YoYInflationIndex>> yoyInflationIndices_;
    map<pair<string, string>, Handle<CPICapFloorTermPriceSurface>> cpiInflationCapFloorPriceSurfaces_;
    map<pair<string, string>, Handle<CPIVolatilitySurface>> cpiInflationCapFloorVolatilitySurfaces_;
    map<pair<string, string>, Handle<YoYCapFloorTermPriceSurface>> yoyInflationCapFloorPriceSurfaces_;
    map<pair<string, string>, Handle<Quote>> equitySpots_;
    map<pair<string, string>, Handle<BlackVolTermStructure>> equityVols_;
    map<pair<string, string>, Handle<Quote>> securitySpreads_;
    map<pair<string, string>, Handle<QuantExt::InflationIndexObserver>> baseCpis_;
    map<tuple<string, string, string>, Handle<QuantExt::CorrelationTermStructure>> correlationCurves_;
    map<pair<string, string>, QuantLib::Handle<QuantLib::Quote>> commoditySpots_;
    map<pair<string, string>, QuantLib::Handle<QuantExt::PriceTermStructure>> commodityCurves_;
    map<pair<string, string>, QuantLib::Handle<QuantLib::BlackVolTermStructure>> commodityVols_;
    map<pair<string, string>, QuantLib::Handle<QuantExt::EquityIndex>> equityCurves_;
    map<pair<string, string>, Handle<Quote>> cprs_;
    Conventions conventions_;

    //! add a swap index to the market
    void addSwapIndex(const string& swapindex, const string& discountIndex,
                      const string& configuration = Market::defaultConfiguration);

    // set of term structure pointers for refresh (per configuration)
    map<string, std::set<boost::shared_ptr<TermStructure>>> refreshTs_;
};
} // namespace data
} // namespace ore