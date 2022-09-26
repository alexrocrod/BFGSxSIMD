#pragma once
#include "boost/math/policies/policy.hpp"
#include <numeric>
#include <boost/math/distributions/normal.hpp>



namespace aadcboost {
    namespace math {


    namespace detail {

        template<class Distribution>
        struct MapAADCBoostToNativeBoost {
            typedef Distribution BOOST_TYPE;
        };


    };

    template <class RealT, class Policy  = ::boost::math::policies::policy<>>
    class chi_squared_distribution: public ::boost::math::chi_squared_distribution<double, Policy> {
    public:
        typedef ::boost::math::chi_squared_distribution<double, Policy> BOOST_BASE;

        chi_squared_distribution(RealT lambda) 
            : ::boost::math::chi_squared_distribution<double, Policy>(AAD_PASSIVE(lambda))
        {}
    };

    namespace detail {
        template<class RealT, class Policy>
        struct MapAADCBoostToNativeBoost<chi_squared_distribution<RealT, Policy>> {
            typedef ::boost::math::chi_squared_distribution<double, Policy> BOOST_TYPE;
        };
    }

    template <class RealT, class Policy = ::boost::math::policies::policy<> >
    class non_central_chi_squared_distribution: public ::boost::math::non_central_chi_squared_distribution<double, Policy> {
    public:
        typedef ::boost::math::non_central_chi_squared_distribution<double, Policy> BOOST_BASE;

        non_central_chi_squared_distribution(RealT df_, RealT lambda) 
            : ::boost::math::non_central_chi_squared_distribution<double, Policy>(AAD_PASSIVE(df_), AAD_PASSIVE(lambda))
        {}
    };

    namespace detail {
        template<class RealT, class Policy>
        struct MapAADCBoostToNativeBoost<non_central_chi_squared_distribution<RealT, Policy>> {
            typedef ::boost::math::non_central_chi_squared_distribution<double, Policy> BOOST_TYPE;
        };
    }


    template <class RealT, class Policy = ::boost::math::policies::policy<> >
    class normal_distribution : public ::boost::math::normal_distribution<double, Policy> {
    public:
        typedef ::boost::math::normal_distribution<double, Policy> BOOST_BASE;
        normal_distribution()
            : ::boost::math::normal_distribution<double, Policy>()
        {}
        normal_distribution(RealT average_, RealT sigma_)
            : ::boost::math::normal_distribution<double, Policy>(AAD_PASSIVE(average_), AAD_PASSIVE(sigma_))
        {}
    };

    namespace detail {
        template<class RealT, class Policy>
        struct MapAADCBoostToNativeBoost<normal_distribution<RealT, Policy>> {
            typedef ::boost::math::normal_distribution<double,Policy> BOOST_TYPE;
        };
    }

    template <class RealT, class Policy = ::boost::math::policies::policy<> >
    class students_t_distribution : public ::boost::math::students_t_distribution<double, Policy> {
    public:
        typedef ::boost::math::students_t_distribution<double, Policy> BOOST_BASE;
        students_t_distribution(RealT df)
            : ::boost::math::students_t_distribution<double, Policy>(AAD_PASSIVE(df))
        {}
    };

    namespace detail {
        template<class RealT, class Policy>
        struct MapAADCBoostToNativeBoost<students_t_distribution<RealT, Policy>> {
            typedef ::boost::math::students_t_distribution<double,Policy> BOOST_TYPE;
        };
    }


template <class Distribution>
inline typename Distribution::value_type pdf(const Distribution& dist, const idouble& x)
{
    typedef typename detail::MapAADCBoostToNativeBoost<Distribution>::BOOST_TYPE boost_dist;
    return ::boost::math::pdf((const boost_dist&)dist, (AAD_PASSIVE(x)));
}

template <class Distribution>
inline typename Distribution::value_type cdf(const Distribution& dist, const idouble& x)
{
   typedef typename detail::MapAADCBoostToNativeBoost<Distribution>::BOOST_TYPE boost_dist;
   return ::boost::math::cdf((const boost_dist&)dist, (AAD_PASSIVE(x)));
}
template <class Distribution>
inline typename Distribution::value_type quantile(const Distribution& dist, const idouble& x)
{
   typedef typename detail::MapAADCBoostToNativeBoost<Distribution>::BOOST_TYPE boost_dist;
   return ::boost::math::quantile((const boost_dist&)dist, (AAD_PASSIVE(x)));
}

inline double lgamma(const idouble& x)
{
   return ::boost::math::lgamma(AAD_PASSIVE(x), 0, ::boost::math::policies::policy<>());
}

inline double
   gamma_p(const idouble& a, const idouble& z)
{
   return ::boost::math::gamma_p(AAD_PASSIVE(a), AAD_PASSIVE(z), ::boost::math::policies::policy<>());
}

inline double
   gamma_p_inv(const idouble& a, const idouble& p)
{
   return ::boost::math::gamma_p_inv(AAD_PASSIVE(a), AAD_PASSIVE(p), ::boost::math::policies::policy<>());
}

inline double
   gamma_q(const idouble& a, const idouble& z)
{
   return ::boost::math::gamma_q(AAD_PASSIVE(a), AAD_PASSIVE(z), ::boost::math::policies::policy<>());
}

inline double
   gamma_q_inv(const idouble& a, const idouble& p)
{
   return ::boost::math::gamma_q_inv(AAD_PASSIVE(a), AAD_PASSIVE(p), ::boost::math::policies::policy<>());
}


} // math
} // aadcboost

namespace aadcboost {
namespace math {

inline ::idouble erf(::idouble z)
{
   return std::erf(z);
}

} // math
} // aadcboost

namespace boost {
    template<>
    class is_arithmetic<idouble> : public is_arithmetic<double> {};
}

