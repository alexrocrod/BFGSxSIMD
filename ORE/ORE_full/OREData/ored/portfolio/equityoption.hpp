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

/*! \file portfolio/equityoption.hpp
    \brief Equity Option data model and serialization
    \ingroup tradedata
*/

#pragma once

#include <ored/portfolio/optiondata.hpp>
#include <ored/portfolio/trade.hpp>

namespace ore {
namespace data {
using std::string;

//! Serializable Equity Option
/*!
  \ingroup tradedata
*/
class EquityOption : public Trade {
public:
    //! Default constructor
    EquityOption() : Trade("EquityOption"), strike_(0.0), quantity_(0.0) {}
    //! Constructor
    EquityOption(Envelope& env, OptionData option, string equityName, string currency, Real strike, Real quantity)
        : Trade("EquityOption", env), option_(option), eqName_(equityName), currency_(currency), strike_(strike),
          quantity_(quantity) {}

    //! Build QuantLib/QuantExt instrument, link pricing engine
    void build(const boost::shared_ptr<EngineFactory>&) override;

    //! Return no fixings for an EquityOption
    std::map<std::string, std::set<QuantLib::Date>> fixings(
        const QuantLib::Date& settlementDate = QuantLib::Date()) const override {
        return {};
    }

    //! \name Inspectors
    //@{
    const OptionData& option() const { return option_; }
    const string& equityName() const { return eqName_; }
    const string& currency() const { return currency_; }
    Real strike() const { return strike_; }
    Real quantity() const { return quantity_; }
    //@}

    //! \name Serialisation
    //@{
    virtual void fromXML(XMLNode* node) override;
    virtual XMLNode* toXML(XMLDocument& doc) override;
    //@}
private:
    OptionData option_;
    string eqName_;
    string currency_;
    Real strike_;
    Real quantity_;
};
} // namespace data
} // namespace ore