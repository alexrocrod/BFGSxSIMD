#pragma once

#include <random>
#include <aadc/idouble.h>

namespace aadcstd {

    template<class RealType>
    class exponential_distribution : public std::exponential_distribution<RealType> {
    };

    template<>
    class exponential_distribution<idouble>  {
    public:
        typedef idouble result_type;
    public:
        exponential_distribution() {}
        explicit exponential_distribution( idouble lambda ) : _u(AAD_PASSIVE(lambda)) {}

        void reset() {
            _u.reset();
        }

        template< class Generator >
        result_type operator()( Generator& g ) {
            return _u(g);
        }

    private:
        std::exponential_distribution<double> _u;
    };

    template<class RealType>
    class gamma_distribution : public std::gamma_distribution<RealType> {
    };

    template<>
    class gamma_distribution<idouble>  {
    public:
        typedef idouble result_type;
    public:
        gamma_distribution() {}
        explicit gamma_distribution( const idouble& alpha, const idouble& beta = 1.0 ) : _u(AAD_PASSIVE(alpha), AAD_PASSIVE(beta)) {}

        void reset() {
            _u.reset();
        }

        template< class Generator >
        result_type operator()( Generator& g ) {
            return _u(g);
        }

    private:
        std::gamma_distribution<double> _u;
    };


    template<class RealType>
    class normal_distribution : public std::normal_distribution<RealType> {
    };

    template<>
    class normal_distribution<idouble>  {
    public:
        typedef idouble result_type;
    public:
        normal_distribution() {}
        explicit normal_distribution( const idouble& mean, const idouble& stddev = 1.0 ) : _u(AAD_PASSIVE(mean), AAD_PASSIVE(stddev)) {}

        void reset() {
            _u.reset();
        }

        template< class Generator >
        result_type operator()( Generator& g ) {
            return _u(g);
        }

    private:
        std::normal_distribution<double> _u;
    };

    template<class RealType>
    class uniform_real_distribution : public std::uniform_real_distribution<RealType> {
    };

    template<>
    class uniform_real_distribution<idouble>  {
    public:
        typedef idouble result_type;
    public:
        uniform_real_distribution() {}
        explicit uniform_real_distribution( const idouble& a, const idouble& b = 1.0 ) : _u(AAD_PASSIVE(a), AAD_PASSIVE(b)) {}

        void reset() {
            _u.reset();
        }

        template< class Generator >
        result_type operator()( Generator& g ) {
            return _u(g);
        }
    private:
        std::uniform_real_distribution<double> _u;
    };


    template<class RealType>
    class cauchy_distribution : public std::cauchy_distribution<RealType> {
    };

    template<>
    class cauchy_distribution<idouble>  {
    public:
        typedef idouble result_type;
    public:
        cauchy_distribution() {}
        explicit cauchy_distribution( const idouble& a, const idouble& b = 1.0 ) : _u(AAD_PASSIVE(a), AAD_PASSIVE(b)) {}

        void reset() {
            _u.reset();
        }

        template< class Generator >
        result_type operator()( Generator& g ) {
            return _u(g);
        }
    private:
        std::cauchy_distribution<double> _u;
    };

}
