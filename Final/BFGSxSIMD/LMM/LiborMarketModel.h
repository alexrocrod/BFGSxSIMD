#ifndef LIBORMARKETMODEL_H
#define LIBORMARKETMODEL_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <aadc/idouble.h>

struct LmmData
{
    std::vector<double>     t;   // diffusion times
    std::vector<double>     T;   // forward times
    double                  a;   // model period in years like 0.25 for 3m Libor model, 0.5 for 6m, etc
    double                  n;   // skew parameter. 1.0 for standard LMM, n > 1 more Hull-Whitish
                                 // simple rate for the period n*a assumed lognormal
};

typedef std::vector<double>::size_type idx_t;


class LmmCalendar
{

    const LmmData &         _lmmData;
    std::vector<idx_t>      _first_fd_num;       // [dd]
    std::vector<idx_t>      _offset;             // [dd]
    idx_t                   _n_dfs;              // number of discount factors per mc path     
    idx_t                   _fd_size;     

    public:
    LmmCalendar(const LmmData & data)
    : _lmmData(data)
    {
        
        _first_fd_num.resize(data.t.size());
        _offset.resize(data.t.size());
        _fd_size = data.T.size();

        _n_dfs = 0;
        for(idx_t dd = 0; dd < data.t.size(); ++dd)
        {
            auto  t_pos = std::lower_bound(data.T.begin(),data.T.end(), data.t[dd]);
            idx_t first_fd = t_pos - data.T.begin();
            _first_fd_num[dd] = first_fd;
            idx_t num_dates = data.T.end() - t_pos; // num fwd dates on this slice
            _offset[dd] = _n_dfs;
            _n_dfs += num_dates;

        }

    }

    idx_t dd_size() const                            { return _first_fd_num.size();}
    idx_t fd_first(idx_t dd) const                   { return _first_fd_num[dd];}
    idx_t fd_size(idx_t dd) const                    { return _fd_size;}
    const LmmData & data() const                     { return _lmmData;}
    idx_t pos(idx_t dd,idx_t fd) const               { return _offset[dd]+fd-_first_fd_num[dd];}
    idx_t num_points() const                         { return _n_dfs;}
    
};

template <class Flt>
class PathData
{
    const LmmCalendar &         _cal;
    std::vector<Flt>            _data;
 
    public:
    Flt & operator () (idx_t dd, idx_t fd)              { return _data[_cal.pos(dd,fd)];}
    const Flt & operator () (idx_t dd, idx_t fd) const  { return _data[_cal.pos(dd,fd)];}
    const LmmCalendar& cal() const                      { return _cal;}

    PathData(const LmmCalendar & calendar)
    : _cal(calendar)
    {
        _data.resize(_cal.num_points(),Flt(0.0));
    }
};

template <class Flt>
class IDiscountCurve
{
    public:
    virtual Flt operator()(double T) const = 0;
};

template <class Flt>
class IVolatilityFunction
{
    public:
    virtual int nFactors() const = 0;
    virtual Flt operator()(double t, double T, int f) const = 0;
};

/*
    Delta in our notation is discounting for single model period between dates in grid T[i].
    In LMM, we model it as a function of the varable x Delta(x) such that

    d x = dW + (some drift)*dt

    $$ Delta(x) = 1 \over {(1 + n a r_0 \exp(x))^(1 \over n) } $$

    what we really need for this implementation, is $ d(ln(\Delta(x))) \over {d x} $, computed at x =\Delta^{-1}(\Delta_0).
    We call it here $$ F(\Delta_0) $$. It is easy to show that 

    \partial {\ln(\Delta)} \over \partial{x} = (-1\over n)(1-\Delta^n)

*/

template <class Flt>
Flt F(Flt Delta0, double n)
{
    return -(1.0 / n)*(1.0-std::pow(Delta0,n));
}

class Rng{
    public:
    virtual double operator() ()= 0;  // N(0,1) RNG sample
};


template <class Flt, class DiscountCurve, class VolatilityFunction>
void propagate(PathData<Flt>& path, const DiscountCurve& discount, 
               const VolatilityFunction& vol, std::vector<Flt>& xi)
{
    const LmmCalendar & cal = path.cal();
    for(idx_t fd = 0; fd < cal.fd_size(0); ++fd)
        path(0,fd) = discount(cal.data().T[fd]);
    
    int n_factors = vol.nFactors();
    double skew_n = cal.data().n;

    std::vector<Flt> zLogVol(n_factors);//, 0.0);          // init with 0.0
    for(idx_t dd = 0; dd < cal.dd_size() - 1; ++dd)
    {
//        if ((dd & 1) && idouble::recording) CAAD_CheckPoint();
        
//        CAAD_CodeFlush();
//        if (idouble::recording) CAAD_LoopPulse(dd);                 // Main loop in time
        std::fill(zLogVol.begin(), zLogVol.end(), 0.);// init with 0.0
        // std::vector<Flt> zLogVol(n_factors, 0.0);          // init with 0.0
        idx_t first_fd = cal.fd_first(dd+1);        
        path(dd+1,first_fd) = path(dd,first_fd);

        double t = cal.data().t[dd];
        double dt = cal.data().t[dd+1]-t;
        double sqrt_dt = std::sqrt(dt);

        for(idx_t fd = first_fd; fd < cal.fd_size(dd) - 1; ++fd)
        {
//            if (idouble::recording) CAAD_CodeFlush();
          if (idouble::recording) CAAD_LoopPulse(fd - first_fd);                 // Main loop in time

            Flt delta0 = path(dd,fd+1)/path(dd,fd);
            Flt ddelta_over_dx = F(delta0,skew_n);
            double T_ = cal.data().T[fd];
            for(idx_t f = 0; f < n_factors; ++f)
            {
                zLogVol[f] += vol(t,T_,f) * ddelta_over_dx ;
            }
            Flt zLogVolSquared = Flt(0.0);
            Flt zLogVol_times_xi = Flt(0.0);
            for(idx_t f = 0; f < n_factors; ++f)
            {
                Flt & v = zLogVol[f];
                zLogVolSquared += v*v;
                zLogVol_times_xi += v * xi[dd*n_factors+f];
            }
            path(dd+1,fd+1) = path(dd,fd+1) * std::exp(zLogVol_times_xi*sqrt_dt-0.5*dt*zLogVolSquared);
        }

    }
}

// pricing of calibration swaptions
// only swaptions with dates on the grid LmmData.T 
// are supported

template <class Flt>
struct cf
{
    idx_t       pmtFd;  // pmt on data.T[pmtFd]
    Flt         pmtAmount;
};

template <class Flt>
struct swaption_data
{
    std::vector<cf<Flt>>    dv01;
    std::vector<cf<Flt>>    flt;
    idx_t               expiryDd; // expiry in data.t[expirtDd] 
};

template <class Flt>
swaption_data<Flt> mock_swaption(double expiry, double tenor, const LmmData & data)
{
    auto t = data.t;
    auto T = data.T;
    auto i_exp = std::lower_bound(t.begin(),t.end(),expiry)-t.begin();
    auto i_beg = std::lower_bound(T.begin(),T.end(), expiry)-T.begin();
    auto i_end = std::lower_bound(T.begin(),T.end(), expiry+tenor)-T.begin();
    swaption_data<Flt> ret;
    ret.flt.resize(2);
    ret.flt[0].pmtFd = i_beg; ret.flt[0].pmtAmount =  1.0;
    ret.flt[1].pmtFd = i_end; ret.flt[1].pmtAmount = -1.0;
    for(idx_t i = i_beg+1; i <= i_end; i++)
    {
        cf<Flt> flow;
        flow.pmtFd=i;
        flow.pmtAmount = data.a;
        ret.dv01.push_back(flow);
    }
    ret.expiryDd = i_exp;
    return ret;
}

template <class Flt>
Flt pv(idx_t asOfDd, const std::vector<cf<Flt>> & flows, const PathData<Flt>& dfs)
{
    Flt ret = 0.0;
    for(idx_t i = 0; i < flows.size(); i++)
    {
        ret += dfs(asOfDd,flows[i].pmtFd)*flows[i].pmtAmount;
    }
    return ret;
}

template <class Flt, class DiscCurve>
Flt pv(const std::vector<cf<Flt>> & flows, const DiscCurve& curve, const LmmData& data)
{
    Flt ret = 0.0;
    for(idx_t i = 0; i < flows.size(); i++)
    {
        ret += curve(data.T[flows[i].pmtFd])*flows[i].pmtAmount;
    }
    return ret;
}

template <class Flt>
std::pair<Flt,Flt> swaption_pv(const swaption_data<Flt>& swaption, const PathData<Flt>& dfs)
{
    idx_t dd = swaption.expiryDd;
    return std::make_pair(pv(dd,swaption.dv01,dfs), pv(dd,swaption.flt,dfs));
}


#endif
