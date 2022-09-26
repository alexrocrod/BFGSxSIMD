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

/*! \file portfolio/fxoption.hpp
    \brief FX Option data model and serialization
    \ingroup tradedata
*/

#pragma once

#include <ored/portfolio/optiondata.hpp>
#include <ored/portfolio/trade.hpp>

namespace ore {
namespace data {
using std::string;

//! Serializable FX Option
/*!
  \ingroup tradedata
*/
class FxOption : public Trade {
public:
    //! Default constructor
    FxOption() : Trade("FxOption"), boughtAmount_(0.0), soldAmount_(0.0) {}
    //! Constructor
    FxOption(Envelope& env, OptionData option, string boughtCurrency, Real boughtAmount, string soldCurrency,
             Real soldAmount)
        : Trade("FxOption", env), option_(option), boughtCurrency_(boughtCurrency), boughtAmount_(boughtAmount),
          soldCurrency_(soldCurrency), soldAmount_(soldAmount) {}

    //! Build QuantLib/QuantExt instrument, link pricing engine
    void build(const boost::shared_ptr<EngineFactory>&) override;

    //! Return no fixings for an FxOption.
    std::map<std::string, std::set<QuantLib::Date>> fixings(
        const QuantLib::Date& settlementDate = QuantLib::Date()) const override {
        return {};
    }

    //! \name Inspectors
    //@{
    const OptionData& option() const { return option_; }
    const string& boughtCurrency() const { return boughtCurrency_; }
    Real boughtAmount() const { return boughtAmount_; }
    const string& soldCurrency() const { return soldCurrency_; }
    Real soldAmount() const { return soldAmount_; }
    //@}

    //! \name Serialisation
    //@{
    virtual void fromXML(XMLNode* node) override;
    virtual XMLNode* toXML(XMLDocument& doc) override;
    //@}
private:
    OptionData option_;
    string boughtCurrency_;
    Real boughtAmount_;
    string soldCurrency_;
    Real soldAmount_;
};
} // namespace data
} // namespace ore
