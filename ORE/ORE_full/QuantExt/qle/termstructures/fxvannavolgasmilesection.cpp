/*
 Copyright (C) 2017 Quaternion Risk Management Ltd
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

#include <ql/math/distributions/normaldistribution.hpp>
#include <qle/termstructures/fxvannavolgasmilesection.hpp>

using namespace QuantLib;

namespace QuantExt {

VannaVolgaSmileSection::VannaVolgaSmileSection(Real spot, Real rd, Real rf, Time t, Volatility atmVol, Volatility rr25d,
                                               Volatility bf25d)
    : FxSmileSection(spot, rd, rf, t, atmVol, rr25d, bf25d) {

    // Consistent Pricing of FX Options
    // Castagna & Mercurio (2006)

    // eq(1). Assumes Delta is unadjusted spot delta which is probably wrong
    k_atm_ = spot * std::exp((rd - rf + 0.5 * atmVol * atmVol) * t);

    // eq(4) + (5).
    vol_25c_ = atmVol + bf25d + 0.5 * rr25d;
    vol_25p_ = atmVol + bf25d - 0.5 * rr25d;

    // eq(6) + (7)
    Real tmp = 0.25 * std::exp(rf * t);
    QL_REQUIRE(tmp > 0 && tmp < 1, "Invalid rate (" << rf << ") and time (" << t << ")in VannaVolgaSmileSection");
    Real alpha = -InverseCumulativeNormal::standard_value(tmp);

    k_25p_ = spot * std::exp((-alpha * vol_25p_ * std::sqrt(t)) + (rd - rf + 0.5 * vol_25p_ * vol_25p_) * t);
    k_25c_ = spot * std::exp((alpha * vol_25c_ * std::sqrt(t)) + (rd - rf + 0.5 * vol_25c_ * vol_25c_) * t);
}

Real VannaVolgaSmileSection::d1(Real x) const {
    return (std::log(spot_ / x) + (rd_ - rf_ + 0.5 * atmVol_ * atmVol_) * t_) / (atmVol_ * std::sqrt(t_));
}

Real VannaVolgaSmileSection::d2(Real x) const {
    return (std::log(spot_ / x) + (rd_ - rf_ - 0.5 * atmVol_ * atmVol_) * t_) / (atmVol_ * std::sqrt(t_));
}

Volatility VannaVolgaSmileSection::volatility(Real k) const {
    QL_REQUIRE(k >= 0, "Non-positive strike (" << k << ")");

    // eq(14). Note sigma = sigma_ATM here.
    Real k1 = k_25p_;
    Real k2 = k_atm_;
    Real k3 = k_25c_;

    // TODO: Cache the (constant) denominator
    Real r1 = std::log(k2 / k) * std::log(k3 / k) / (std::log(k2 / k1) * std::log(k3 / k1));
    Real r2 = std::log(k / k1) * std::log(k3 / k) / (std::log(k2 / k1) * std::log(k3 / k2));
    Real r3 = std::log(k / k1) * std::log(k / k2) / (std::log(k3 / k1) * std::log(k3 / k2));

    Real sigma1_k = r1 * vol_25p_ + r2 * atmVol_ + r3 * vol_25c_;

    Real D1 = sigma1_k - atmVol_;

    // No middle term as sigma = sigma_atm
    Real D2 = r1 * d1(k1) * d2(k1) * (vol_25p_ - atmVol_) * (vol_25p_ - atmVol_) +
              r3 * d1(k3) * d2(k3) * (vol_25c_ - atmVol_) * (vol_25c_ - atmVol_);

    Real d1d2k = d1(k) * d2(k);

    Real tmp = atmVol_ * atmVol_ + d1d2k * (2 * atmVol_ * D1 + D2);
    QL_REQUIRE(tmp >= 0, "VannaVolga attempting to take square root of negative number");

    return atmVol_ + (-atmVol_ + std::sqrt(tmp)) / d1d2k;
}

} // namespace QuantExt
