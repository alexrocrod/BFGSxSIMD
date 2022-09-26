#ifndef LMM_TEST_DATA_H
#define LMM_TEST_DATA_H

// mock data to test implementation

#include "LiborMarketModel.h"

LmmData test_data();

template <class Flt>
class TestCurve // : public IDiscountCurve<Flt>
{
    public:
    Flt  operator()(double T) const
    {
        return std::exp( - T * 0.05); // 5% flat
    }
};

template <class Flt>
class TestVol: public IVolatilityFunction<Flt>
{
    public:
    virtual int nFactors() const            { return 1;}
    virtual Flt operator()(double t, double T, int f) const 
    {
        return 0.05 + 0.1 * std::exp(-(T-t)/5.0);
    }
};


// simple 3 params
template <class Flt>
class LmmVol1
{
    public:

    std::vector<Flt> p;

    LmmVol1()
    {
        p.resize(nParameters());
        setStartingPoint0();
    }

    void setStartingPoint0()
    {
        // some starting point to compute 
        // target prices for calibration
        p[0] = 0.1;   // 10% long term vol
        p[1] = 1.0;   // 100% increase for short maturities
        p[2] = 0.2;   // 5y time scale for increase
    }

    void setStartingPoint1()
    {
        // starting point to start test calibratuion
        setStartingPoint0();
        for(int i = 0; i < p.size(); ++i)
            p[i] *= 0.9;
    }

    
    static int nFactors()                { return 1;}
    static int nParameters()             { return 3;}

    static std::vector<double> const & params_from()
    {
        static std::vector<double> ret(nParameters());
        ret[0] = 0.001;
        ret[1] = 0.05;
        ret[2] = 0.02;
        return ret;
    }

    static std::vector<double> const & params_to()
    {
        static std::vector<double> ret(nParameters());
        ret[0] = 2.0;
        ret[1] = 5.0;
        ret[2] = 5.0;
        return ret;
    }
    
    Flt operator()(double t, double T, int f) const 
    {
        return p[0] * (1.0 + p[1] * std::exp(- p[2] * (T-t)));
    }
};

// Rebonato
template <class Flt>
class LmmVolR
{
        public:

    std::vector<Flt> p;

    LmmVolR()
    {
        p.resize(nParameters());
        setStartingPoint0();
    }

    void setStartingPoint0()
    {
        // some starting point to compute 
        // target prices for calibration
        p[0] = 0.1;   // 10% long term vol
        p[1] = 1.0;   // 100% increase for short maturities
        p[2] = 0.2;
        p[3] = 0.2;   // 5y time scale for increase
        p[4] = 1.0;   // 100% increase for short maturities
        p[5] = 0.2;
        p[6] = 0.2;   // 5y time scale for increase
    }

    void setStartingPoint1()
    {
        // starting point to start test calibratuion
        setStartingPoint0();
        for(int i = 0; i < p.size(); ++i)
            p[i] *= 0.9;
    }

    
    static int nFactors()                { return 1;}
    static int nParameters()             { return 7;}

    static std::vector<double> const & params_from()
    {
        static std::vector<double> ret(nParameters());
        ret[0] = 0.001;
        ret[1] = 0.05;
        ret[2] = 0.0;
        ret[3] = 0.02;
        ret[4] = 0.05;
        ret[5] = 0.0;
        ret[6] = 0.02;
        return ret;
    }

    static std::vector<double> const & params_to()
    {
        static std::vector<double> ret(nParameters());
        ret[0] = 2.0;
        ret[1] = 5.0;
        ret[2] = 5.0;
        ret[3] = 5.0;
        ret[4] = 5.0;
        ret[5] = 5.0;
        ret[6] = 5.0;
        return ret;
    }
    
    Flt operator()(double t, double T, int f) const 
    {
        Flt vol =  p[0] * (1.0 + (p[1] + p[2] * (T-t)) * std::exp(- p[3] * (T-t))) *
                      (1.0 + (p[4] + p[5] * (T  )) * std::exp(- p[6] * (T  )));
        return vol;
    }

};

// select which one to use:
template <class Flt>
using LmmVol = LmmVolR<Flt>;


std::vector<std::pair<double,double>>   test_swaptions();

std::vector<std::pair<double,double>>   std_swaptions(const LmmData& data);


#endif
