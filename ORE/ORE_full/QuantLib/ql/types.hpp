/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005 StatPro Italia srl

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file types.hpp
    \brief Custom types
*/

#ifndef quantlib_types_hpp
#define quantlib_types_hpp

// #define QL_AAD

// #define AADC_QL

#include <cstdint>
#define QL_AAD_REPORT_ERROR_USE throw "Not supported AAD code called";

#ifdef AADC_QL
#include <aadc/idouble.h>
#include <aadc/ibool.h>
#include <aadc/iint.h>
// #include <aadc/aadc_debug.h>
#else
#include <aadc/aadc_compat.h>
#endif

#include <ql/qldefines.hpp>
#include <cstddef>

#include <numeric>
#include <boost/math/distributions/normal.hpp>

namespace QuantLib {

    //! integer number
    /*! \ingroup types */
    typedef QL_INTEGER Integer;

    //! large integer number
    /*! \ingroup types */
    typedef QL_BIG_INTEGER BigInteger;

    //! positive integer
    /*! \ingroup types */
    typedef unsigned QL_INTEGER Natural;

    //! large positive integer
    typedef unsigned QL_BIG_INTEGER BigNatural;

    //! real number
    /*! \ingroup types */
//    typedef QL_REAL Real;
#ifdef AADC_QL
    typedef idouble Real;
#else
    typedef double Real;
#endif 

#ifdef AADC_QL
    typedef ibool QLBool;
#else
    typedef bool QLBool;
#endif 

    //! decimal number
    /*! \ingroup types */
    typedef Real Decimal;

    //! size of a container
    /*! \ingroup types */
    typedef std::size_t Size;

    //! continuous quantity with 1-year units
    /*! \ingroup types */
    typedef Real Time;

    //! discount factor between dates
    /*! \ingroup types */
    typedef Real DiscountFactor;

    //! interest rates
    /*! \ingroup types */
    typedef Real Rate;

    //! spreads on interest rates
    /*! \ingroup types */
    typedef Real Spread;

    //! volatility
    /*! \ingroup types */
    typedef Real Volatility;

    //! probability
    /*! \ingroup types */
    typedef Real Probability;

}

#ifndef AADC_QL
#define AAD_PASSIVE(r) (r)
#endif

namespace QuantLib {
    inline Integer Integer2(const Real& r) {
        return AAD_PASSIVE(r);
    }
    inline std::size_t Size2(const Real& r) {
        return AAD_PASSIVE(r);
    }
    inline BigInteger BigInteger2(const Real& r) {
        return AAD_PASSIVE(r);
    }
    inline BigNatural BigNatural2(const Real& r) {
        return AAD_PASSIVE(r);
    }
    
}

#ifdef AADC_QL
namespace std {
  template<typename _InputIterator1, typename _InputIterator2>
    inline QuantLib::Real
    inner_product(_InputIterator1 __first1, _InputIterator1 __last1,
		  _InputIterator2 __first2, double __init)
    {
      return inner_product(__first1, __last1, __first2, QuantLib::Real(__init));
    }


  template<typename _InputIterator>
    inline QuantLib::Real
    accumulate(_InputIterator __first, _InputIterator __last, double __init)
    {
      return accumulate(__first, __last, QuantLib::Real(__init));
    }

}
#endif // QL_AAD

namespace stdpassive {
#ifdef AADC_QL
    inline double floor(const idouble& val) {
        return std::floor(AAD_PASSIVE(val));
    }
    inline double ceil(const idouble& val) {
        return std::ceil(AAD_PASSIVE(val));
    }
#endif
    inline double floor(double val) {
        return std::floor(val);
    }
    inline double ceil(double val) {
        return std::ceil(val);
    }
}

#ifdef AADC_QL

#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/type_traits/is_floating_point.hpp>

namespace boost {
   template<> struct is_floating_point<idouble> : public true_type{};

   template<> 
   inline bool boost::math::isfinite<idouble>(idouble a) { return boost::math::isfinite<double>(AAD_PASSIVE(a)); }
};

namespace aadc {
namespace boost {
    namespace math {


    namespace detail {

        template<class Distribution>
        struct MapAADCBoostToNativeBoost {
            typedef Distribution BOOST_TYPE;
        };


    };

    template <class RealT>
    class chi_squared_distribution: public ::boost::math::chi_squared_distribution<double> {
    public:
        typedef ::boost::math::chi_squared_distribution<double> BOOST_BASE;

        chi_squared_distribution(RealT lambda) 
            : ::boost::math::chi_squared_distribution<double>(AAD_PASSIVE(lambda))
        {}
    };

    namespace detail {
        template<class RealT>
        struct MapAADCBoostToNativeBoost<chi_squared_distribution<RealT>> {
            typedef ::boost::math::chi_squared_distribution<double> BOOST_TYPE;
        };
    }

    template <class RealT>
    class non_central_chi_squared_distribution: public ::boost::math::non_central_chi_squared_distribution<double> {
    public:
        typedef ::boost::math::non_central_chi_squared_distribution<double> BOOST_BASE;

        non_central_chi_squared_distribution(RealT df_, RealT lambda) 
            : ::boost::math::non_central_chi_squared_distribution<double>(AAD_PASSIVE(df_), AAD_PASSIVE(lambda))
        {}
    };

    namespace detail {
        template<class RealT>
        struct MapAADCBoostToNativeBoost<non_central_chi_squared_distribution<RealT>> {
            typedef ::boost::math::non_central_chi_squared_distribution<double> BOOST_TYPE;
        };
    }


    template <class RealT>
    class normal_distribution : public ::boost::math::normal_distribution<double> {
    public:
        typedef ::boost::math::normal_distribution<double> BOOST_BASE;
        normal_distribution(RealT average_, RealT sigma_)
            : ::boost::math::normal_distribution<double>(AAD_PASSIVE(average_), AAD_PASSIVE(sigma_))
        {}
    };

    namespace detail {
        template<class RealT>
        struct MapAADCBoostToNativeBoost<normal_distribution<RealT>> {
            typedef ::boost::math::normal_distribution<double> BOOST_TYPE;
        };
    }

template <class Distribution>
inline typename Distribution::value_type pdf(const Distribution& dist, const QuantLib::Real& x)
{
    typedef typename detail::MapAADCBoostToNativeBoost<Distribution>::BOOST_TYPE boost_dist;
    return ::boost::math::pdf((const boost_dist&)dist, (AAD_PASSIVE(x)));
}

template <class Distribution>
inline typename Distribution::value_type cdf(const Distribution& dist, const QuantLib::Real& x)
{
   typedef typename detail::MapAADCBoostToNativeBoost<Distribution>::BOOST_TYPE boost_dist;
   return ::boost::math::cdf((const boost_dist&)dist, (AAD_PASSIVE(x)));
}
template <class Distribution>
inline typename Distribution::value_type quantile(const Distribution& dist, const QuantLib::Real& x)
{
   typedef typename detail::MapAADCBoostToNativeBoost<Distribution>::BOOST_TYPE boost_dist;
   return ::boost::math::quantile((const boost_dist&)dist, (AAD_PASSIVE(x)));
}

inline double lgamma(const QuantLib::Real& x)
{
   return ::boost::math::lgamma(AAD_PASSIVE(x), 0, ::boost::math::policies::policy<>());
}

inline double
   gamma_p(const QuantLib::Real& a, const QuantLib::Real& z)
{
   return ::boost::math::gamma_p(AAD_PASSIVE(a), AAD_PASSIVE(z), ::boost::math::policies::policy<>());
}

inline double
   gamma_p_inv(const QuantLib::Real& a, const QuantLib::Real& p)
{
   return ::boost::math::gamma_p_inv(AAD_PASSIVE(a), AAD_PASSIVE(p), ::boost::math::policies::policy<>());
}

inline double
   gamma_q(const QuantLib::Real& a, const QuantLib::Real& z)
{
   return ::boost::math::gamma_q(AAD_PASSIVE(a), AAD_PASSIVE(z), ::boost::math::policies::policy<>());
}

inline double
   gamma_q_inv(const QuantLib::Real& a, const QuantLib::Real& p)
{
   return ::boost::math::gamma_q_inv(AAD_PASSIVE(a), AAD_PASSIVE(p), ::boost::math::policies::policy<>());
}


} // math
} // boost
} // aadc
#define AADBOOST aadc::boost

#else
#define AADBOOST boost

#endif // QL_AAD

#ifdef AADC_QL
namespace boost {
namespace math {

inline ::QuantLib::Real erf(::QuantLib::Real z)
{
   return std::erf(z);
}

} // math
} // boost
#endif // QL_AAD


// AADC patch. Don't differentiate boost distributions for now
namespace QuantLib {
typedef double BoostReal;
};


#define BOOST_CHECK_CLOSE_PASSIVE( L, R, T ) BOOST_CHECK_CLOSE( AAD_PASSIVE(L), AAD_PASSIVE(R), AAD_PASSIVE(T) ) 
#define BOOST_CHECK_SMALL_PASSIVE( FPV, T ) BOOST_CHECK_SMALL( AAD_PASSIVE(FPV), AAD_PASSIVE(T) )


#endif
