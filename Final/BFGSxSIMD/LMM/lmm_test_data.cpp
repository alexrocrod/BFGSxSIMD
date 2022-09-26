#include "LiborMarketModel.h"
#include "lmm_test_data.h"
#include <set>


LmmData test_data()
{
    LmmData ret;
    ret.a = 0.25;
    ret.n = 1.0;
    ret.T.resize(41);  // 10 years 41 fwd dates
    for(idx_t i = 0; i < ret.T.size(); ++i)
    {
        ret.T[i] = ret.a * i;
    }

    ret.t.reserve(25);
    ret.t.resize(21);  // 5 years 20 3m steps 
    for(idx_t i = 0; i < ret.t.size(); ++i)
    {
        ret.t[i] = ret.a * i;
    }

    ret.t.push_back(ret.a*0.5);
    ret.t.push_back(ret.a*1.5);
    std::sort(ret.t.begin(),ret.t.end()); 

    return ret;

}

std::vector<std::pair<double,double>>
test_swaptions()
{
    
    std::vector<std::pair<double,double>> ret;

    ret.push_back(std::make_pair(1.0,1.0));         // 1 x 1
    ret.push_back(std::make_pair(2.0,1.0));         // 2 x 1  
    ret.push_back(std::make_pair(3.0,1.0));
    ret.push_back(std::make_pair(4.0,1.0));

    return ret;
}


std::vector<std::pair<double,double>>
std_swaptions(const LmmData& data)
{
    
    double std_tenors[] = { 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 15.0, 20.0, 30.0 };
    int n_tenors = sizeof(std_tenors)/sizeof(double);

    std::vector<std::pair<double,double>> ret;
    std::set<double> t_(data.t.begin(),data.t.end());
    std::set<double> T_(data.T.begin(),data.T.end());

    for(int i = 0; i < n_tenors; i++)
    {
        for(int j = 0; j < n_tenors; j++)
        {
            double T0 = std_tenors[i];
            double T1 = std_tenors[j];
            if( (t_.find(T0) != t_.end()) && (T_.find(T0) != T_.end()) && (T_.find(T0+T1)) != T_.end())
                ret.push_back(std::make_pair(T0,T1));
            
        }
    }
    
    return ret;
}

