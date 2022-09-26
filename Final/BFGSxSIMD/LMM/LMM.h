#if (!defined(LMM_H) && !defined(AAD_PASS)) || (!defined(LMM_H_AAD) && defined(AAD_PASS))
#ifndef AAD_PASS
#define LMM_H
#else
#define LMM_H_AAD
#endif

#include <vector>

#include <aadc/idouble.h>
#include <math.h>


class Surface {
public:
     Surface(const std::vector<double>& x, const std::vector<double>& y)
        : m_x(x), m_y(y)
    {}

    void initToConst(const mdouble& val) {
        m_vals.resize(m_x.size());
        for (int xi = 0; xi < m_x.size(); ++xi) {
            m_vals[xi].resize(m_y.size());
            std::fill(m_vals[xi].begin(), m_vals[xi].end(), val);
        }
    }

    mdouble operator ()(const double& x, const double& y) {
        // todo
        return m_vals[0][0];
    }

    std::vector<std::vector<mdouble> > m_vals;
    const std::vector<double> m_x, m_y;
};


class LMM {
public:
    LMM(
        const std::vector<double>& tau,
        const Surface& vol, const mdouble& beta
    ) : m_tau(tau)
        , m_vol(vol)
        , m_beta(beta)
    {}

    template<class Curve>
    void initT0Fwds(const Curve& crv) {
        m_fwds.resize(m_tau.size() - 1);
        for (int i = 0; i < m_fwds.size(); ++i) {
            m_fwds[i] = (crv.df(m_tau[i]) / crv.df(m_tau[i+1]) - 1.0) / (m_tau[i+1] - m_tau[i]);
            m_fwds_init[i] = m_fwds[i];
        }
        m_stoch_vol = 1.0;
    }

    void nextT(const std::vector<mdouble>& normals, const double& dt) {
        int cur_fwd_indx = fwdIndex(m_cur_t + dt);
        
        std::vector<mdouble> drift(normals.size(), 0.0);
        
        for (int fi = cur_fwd_indx; fi < m_fwds.size(); ++fi) {
            mdouble vol = std::sqrt(m_vol(m_cur_t, m_tau[fi] - m_cur_t));
            mdouble sqrt_stoch_vol = std::sqrt(m_stoch_vol);

            mdouble stoch_sum = 0;
            std::vector<mdouble> factor(normals.size(), 0.0);
            for (int j = 0; j < normals.size(); ++j) {
                stoch_sum = stoch_sum + vol * factor[j] * (drift[j] * sqrt_stoch_vol * dt + normals[0] * std::sqrt(dt));
            }
            mdouble fwd_d = (m_beta * m_fwds[fi] + (1 - m_beta) * m_fwds_init[fi]) * sqrt_stoch_vol * stoch_sum;

            // drift = drift + (tau) * vol * m_fwds[fi] / (1 + tau * m_fwds[fi])
            m_fwds[fi] = m_fwds[fi] + fwd_d;
        }
    }
    
    int fwdIndex(const double& t) {
        int res(0);
        while(res < m_tau.size() && (m_tau[res] < t)) ++res;
        return res;
    }

private:
    double m_cur_t;
    mdouble m_stoch_vol;
    std::vector<mdouble> m_fwds;
    std::vector<mdouble> m_fwds_init;
    std::vector<double> m_tau;
    Surface m_vol;
    mdouble m_beta;

};



#ifndef AAD_PASS 
namespace AAD {
    #define AAD_PASS
    #include "LMM.h"
};
#endif

#endif
