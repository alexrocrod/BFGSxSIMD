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

/*! \file portfolio/capfloor.hpp
    \brief Ibor cap, floor or collar trade data model and serialization
    \ingroup tradedata
*/

#pragma once

#include <ored/portfolio/legdata.hpp>
#include <ored/portfolio/trade.hpp>

namespace ore {
namespace data {

//! Serializable cap, floor, collar
/*! \ingroup tradedata
 */
class CapFloor : public Trade {
public:
    CapFloor() : Trade("CapFloor") {}
    CapFloor(const Envelope& env, const string& longShort, const LegData& leg, const vector<Real>& caps,
             const vector<Real>& floors)
        : Trade("CapFloor", env), longShort_(longShort), legData_(leg), caps_(caps), floors_(floors) {}

    virtual void build(const boost::shared_ptr<EngineFactory>&) override;
    
    //! Return the fixings that will be requested to price the CapFloor given the \p settlementDate.
    std::map<std::string, std::set<QuantLib::Date>> fixings(
        const QuantLib::Date& settlementDate = QuantLib::Date()) const override;

    //! Inspectors
    //@{
    const string& longShort() const { return longShort_; }
    const LegData& leg() const { return legData_; }
    const vector<Real>& caps() const { return caps_; }
    const vector<Real>& floors() const { return floors_; }
    //@}

    virtual void fromXML(XMLNode* node) override;
    virtual XMLNode* toXML(XMLDocument& doc) override;

private:
    string longShort_;
    LegData legData_;
    vector<Real> caps_;
    vector<Real> floors_;

    //! Store the name of the underlying floating leg index
    std::string underlyingIndex_;
};
} // namespace data
} // namespace ore