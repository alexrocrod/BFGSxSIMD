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

/*! \file portfolio/fxswap.hpp
    \brief FX Swap data model and serialization
    \ingroup tradedata
*/

#pragma once

#include <ored/portfolio/trade.hpp>

namespace ore {
namespace data {
using std::string;

//! Serializable FX Swap
/*!
  \ingroup tradedata
*/
class FxSwap : public Trade {
public:
    //! Default constructor
    FxSwap()
        : Trade("FxSwap"), nearBoughtAmount_(0.0), nearSoldAmount_(0.0), farBoughtAmount_(0.0), farSoldAmount_(0.0) {}
    //! Constructor
    FxSwap(Envelope& env, string nearDate, string farDate, string nearBoughtCurrency, Real nearBoughtAmount,
           string nearSoldCurrency, Real nearSoldAmount, Real farBoughtAmount, Real farSoldAmount)
        : Trade("FxSwap", env), nearDate_(nearDate), farDate_(farDate), nearBoughtCurrency_(nearBoughtCurrency),
          nearBoughtAmount_(nearBoughtAmount), nearSoldCurrency_(nearSoldCurrency), nearSoldAmount_(nearSoldAmount),
          farBoughtAmount_(farBoughtAmount), farSoldAmount_(farSoldAmount) {}
    /*! Constructs a composite pricing engine of two FX forward pricing engines.
    One with the near amounts as notionals, the other with the far amounts.
    NPV is the total npv of these trades.
    */

    //! Build QuantLib/QuantExt instrument, link pricing engine
    void build(const boost::shared_ptr<EngineFactory>&) override;

    //! Return no fixings for an FxSwap.
    std::map<std::string, std::set<QuantLib::Date>> fixings(
        const QuantLib::Date& settlementDate = QuantLib::Date()) const override {
        return {};
    }

    //! \name Inspectors
    //@{
    string nearDate() { return nearDate_; }
    string farDate() { return farDate_; }
    string nearBoughtCurrency() { return nearBoughtCurrency_; }
    Real nearBoughtAmount() { return nearBoughtAmount_; }
    string nearSoldCurrency() { return nearSoldCurrency_; }
    Real nearSoldAmount() { return nearSoldAmount_; }
    Real farBoughtAmount() { return farBoughtAmount_; }
    Real farSoldAmount() { return farSoldAmount_; }
    //@}

    //! \name Serialisation
    //@{
    virtual void fromXML(XMLNode* node) override;
    virtual XMLNode* toXML(XMLDocument& doc) override;
    //@}
private:
    string nearDate_;
    string farDate_;
    string nearBoughtCurrency_; // farBoughtCurrency==nearSoldCurrency
    Real nearBoughtAmount_;
    string nearSoldCurrency_;
    Real nearSoldAmount_;
    Real farBoughtAmount_;
    Real farSoldAmount_;
};
} // namespace data
} // namespace ore
