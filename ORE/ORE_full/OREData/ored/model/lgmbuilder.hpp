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

/*! \file ored/model/lgmbuilder.hpp
    \brief Build an lgm model
    \ingroup models
*/

#pragma once

#include <map>
#include <ostream>
#include <vector>

#include <qle/models/lgm.hpp>

#include <ored/model/irlgmdata.hpp>
#include <ored/model/modelbuilder.hpp>

namespace ore {
namespace data {
using namespace QuantLib;

//! Builder for a Linear Gauss Markov model component
/*!
  This class is a utility that turns a Linear Gauss Markov
  model description into an interest rate model parametrisation which
  can be used to instantiate a CrossAssetModel.

  \ingroup models
 */
class LgmBuilder : public ModelBuilder {
public:
    /*! The configuration should refer to the calibration configuration here,
      alternative discounting curves are then usually set in the pricing
      engines for swaptions etc. */
    //LgmBuilder(const boost::shared_ptr<ore::data::Market>& market, const boost::shared_ptr<IrLgmData>& data,
    //           const std::string& configuration = Market::defaultConfiguration, Real bootstrapTolerance = 0.001);


    LgmBuilder(int& aadc_i,   std::map<std::tuple<string, YieldCurveType, std::string>, AADCCurveData>& curves_data,
        const boost::shared_ptr<ore::data::Market>& market, const boost::shared_ptr<IrLgmData>& data,
               const std::string& configuration = Market::defaultConfiguration, Real bootstrapTolerance = 0.001
               , bool use_aadc=false
               , bool shift_=false
               );

    //! Return calibration error
    Real error() {
        calculate();
        return error_;
    }

    //! \name Inspectors
    //@{
    std::string currency() { return data_->ccy(); }
    boost::shared_ptr<QuantExt::LGM>& model() {
        calculate();
        return model_;
    }
    boost::shared_ptr<QuantExt::IrLgm1fParametrization>& parametrization() { return parametrization_; }
    RelinkableHandle<YieldTermStructure> discountCurve() { return discountCurve_; }
    std::vector<boost::shared_ptr<BlackCalibrationHelper>> swaptionBasket() {
        calculate();
        return swaptionBasket_;
    }
    //@}
        mutable bool aadc_initialized_;
        bool use_aadc_;
        int& aadc_i_;
        mutable bool shift_=false;      
                std::map<std::tuple<string, YieldCurveType, std::string>, AADCCurveData>& curves_data_;
private:
    void performCalculations() const override;
    void buildSwaptionBasket() const;

    boost::shared_ptr<ore::data::Market> market_;
    const std::string configuration_;
    boost::shared_ptr<IrLgmData> data_;
    Real bootstrapTolerance_;
    mutable Real error_;
    boost::shared_ptr<QuantExt::LGM> model_;
    Array params_;
    boost::shared_ptr<QuantLib::PricingEngine> swaptionEngine_;
    boost::shared_ptr<QuantExt::IrLgm1fParametrization> parametrization_;
    RelinkableHandle<YieldTermStructure> discountCurve_;
    mutable std::vector<boost::shared_ptr<BlackCalibrationHelper>> swaptionBasket_;
    mutable Array swaptionExpiries_;
    mutable Array swaptionMaturities_;

    Handle<QuantLib::SwaptionVolatilityStructure> svts_;
    Handle<SwapIndex> swapIndex_, shortSwapIndex_;

    // TODO: Move CalibrationErrorType, optimizer and end criteria parameters to data
    boost::shared_ptr<OptimizationMethod> optimizationMethod_;
    EndCriteria endCriteria_;
    BlackCalibrationHelper::CalibrationErrorType calibrationErrorType_;
};
} // namespace data
} // namespace ore
