#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#include <set>
#include "LiborMarketModel.h"

#include <aadc/aadc.h>
// #define AADC_INTERNAL_DEBUG
#ifdef AADC_INTERNAL_DEBUG
#include <aadc_private.h>
#endif
#include <aadc/ibool.h>

#include "lmm_test_data.h"
#include <thread>

#include <Eigen/Core>
#include <LBFGS.h>

using Eigen::VectorXd;
using namespace LBFGSpp;
using namespace aadc;

class StdRng: public Rng
{
    std::mt19937                _gen;
    std::normal_distribution<>  _dist;
    public:
    StdRng():
        _dist(0,1)
    {
        // same seed on each run
        _gen.seed(7777);
    }
    virtual double operator() ()
    {
        double ret = _dist(_gen);
        return ret;
    }
};


template <class Flt, class DiscCurve, class LmmVol>
class calibration_helper
{
    public:
    const LmmData&              data;
    LmmCalendar                 cal;
    LmmVol                      vol;
    const DiscCurve&  curve;
    PathData<Flt>               path;

    std::vector<std::pair<double,double>> calibration_swaptions_info;
    std::vector<swaption_data<Flt>>     calibration_swaptions;
    std::vector<Flt>                    calibration_swaptions_rate;  // off the curve
    std::vector<Flt>                    calibration_swaptions_dv01;  // off the curve

    std::vector<Flt>                    path_flt;  // calib swaption floating leg on the path
    std::vector<Flt>                    path_fxd;  // calib swaption fixed leg on the path

    // ---------------- these vaiables are to be tweakable inputs
    std::vector<Flt>        xi;             // N(0,1) for single path, nFactors * nSteps

    // ----------------


    // ------------------------ result of the sinle path run that is needed for calibration
    std::vector<Flt>                    path_calib_swaption_payoff;
    // ------------------------


    calibration_helper(const LmmData& lmmData, 
    const DiscCurve& crv) 
    : data(lmmData), cal(lmmData),curve(crv),path(cal)
    {

        calibration_swaptions_info = std_swaptions(data);
        idx_t nCalibSwaptions = calibration_swaptions_info.size();

        calibration_swaptions.resize(nCalibSwaptions);
        calibration_swaptions_rate.resize(nCalibSwaptions);
        calibration_swaptions_dv01.resize(nCalibSwaptions);

        path_flt.resize(nCalibSwaptions);
        path_fxd.resize(nCalibSwaptions);
        path_calib_swaption_payoff.resize(nCalibSwaptions);

        for(idx_t k = 0; k < calibration_swaptions.size(); k++)
        {
            auto & cs_info = calibration_swaptions_info[k];
            calibration_swaptions[k] = mock_swaption<Flt>(cs_info.first, cs_info.second, data);
            auto & cs = calibration_swaptions[k];
            Flt flt0 = pv(cs.flt,curve,data);
            Flt fxd0 = pv(cs.dv01,curve,data);
            calibration_swaptions_rate[k] = flt0/fxd0;
            calibration_swaptions_dv01[k] = fxd0;
        }

    }


    void run_one_path()
    {

        auto nFactors = vol.nFactors();
        propagate(path,curve,vol,xi);
//        if (idouble::recording) CAAD_CheckPoint();

        idx_t nCalibSwaptions = calibration_swaptions.size();

        for(idx_t k = 0; k < nCalibSwaptions; k++)
        {
//            if (idouble::recording)  CAAD_LoopPulse(k);
            auto & swaption_data = calibration_swaptions[k];
            path_fxd[k] = pv(swaption_data.expiryDd,swaption_data.dv01,path);
            path_flt[k] = pv(swaption_data.expiryDd,swaption_data.flt,path);
        }
//        if (idouble::recording) CAAD_CheckPoint();

        for(idx_t k = 0; k < nCalibSwaptions; k++)
        {     
//            if (idouble::recording) CAAD_LoopPulse(k);
            Flt atm_rate = calibration_swaptions_rate[k];
            Flt swap_pv = path_flt[k] - atm_rate * path_fxd[k];
            
            // here, calibrate to ATM RTP swaptions. collect pv for these
            // note: requires ibool.h
            path_calib_swaption_payoff[k] = std::max(swap_pv,0.0);
            //path_calib_swaption_payoff[k] = swap_pv >= Flt(0.0) ? swap_pv : Flt(0.0);
        }

    } // run_one_path

};

template<class mmType>
class calibration_evaluator
{

    LmmData&    data_;
    TestCurve<idouble> & disc_;
    int nSteps;
    int nPaths;
    int nFactors;
    mmVector<mmType> mmRandomSample;

    std::vector<aadc::AADCArgument> p_idx;
    std::vector<aadc::AADCArgument> xi_idx;
    std::chrono::microseconds simulate1_time;

    int num_ci;
    std::vector<aadc::AADCResult> path_calib_swaption_payoff_idx;

#ifdef AADC_INTERNAL_DEBUG
    std::shared_ptr<aadc::AADCCompiledFunctions<mmType> >   aad_funcs; 
#else
    std::shared_ptr<aadc::AADCFunctions<mmType> >   aad_funcs; 
#endif 
    std::vector<double> ci_pv;
    int nSIMD;
    int nRepeats;
    int nThreads;

    std::vector<double> target_pv;
    std::vector<double> grad;
    double factor = 1000000.0;

    // --added-- 
    bool onlyF = true; // enables the use of only 1D derivatives
    int ncalls_forward = 0;
    int ncalls_reverse = 0;
    double delta = 1e-7; // default: 1e-7
    // ---------

    public:

    // Constructor
    calibration_evaluator(
        LmmData& data, TestCurve<idouble> & disc, TestCurve<double> & disc_org, int num_threads)
        : data_(data),disc_(disc)
    {
        LmmCalendar         cal(data_);
        StdRng              gen;

        nSteps = data.t.size() - 1;
        nPaths = 1024*64;
        nFactors = 1;
        nThreads = num_threads;
        nSIMD = sizeof(mmType) / sizeof(double);

        mmRandomSample.resize(nSteps*nFactors*nPaths / nSIMD);
        double* randomSample = (double*)(&(mmRandomSample[0]));
        {
            auto gen_start = std::chrono::high_resolution_clock::now(); 
            for(idx_t p = 0; p < nSteps*nFactors*nPaths; p++)
                randomSample[p] = gen();
            auto gen_stop = std::chrono::high_resolution_clock::now(); 
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(gen_stop - gen_start); 
            std::cout <<" generating randoms :" << duration.count() << std::endl;
        }
        calibration_helper<idouble, TestCurve<idouble>, LmmVol<idouble>> lmm_calibration_helper(data, disc);

        auto compile_start = std::chrono::high_resolution_clock::now();

#ifdef AADC_INTERNAL_DEBUG
        aad_funcs = std::shared_ptr<aadc::AADCCompiledFunctions<mmType> >(
            new aadc::AADCCompiledFunctions<mmType>(
                16|8|4, "examples/LMM/functions/", "LmmOnePath"
            )
        );
#else
        aad_funcs = std::shared_ptr<aadc::AADCFunctions<mmType> >(
            new aadc::AADCFunctions<mmType>()
        );
#endif 

        lmm_calibration_helper.xi.resize(nSteps*nFactors);


        num_ci = lmm_calibration_helper.calibration_swaptions_rate.size();
        path_calib_swaption_payoff_idx.resize(num_ci);


        aad_funcs->startRecording();
        p_idx.resize(lmm_calibration_helper.vol.nParameters());
        for(int i = 0; i < p_idx.size(); i++)
            p_idx[i] = lmm_calibration_helper.vol.p[i].markAsInput();
        
        xi_idx.resize(nSteps*nFactors);
        for(idx_t i = 0; i < xi_idx.size(); i++)
            xi_idx[i]=lmm_calibration_helper.xi[i].markAsInputNoDiff();
        
        
        for (int k = 0; k < 1; ++k)
        lmm_calibration_helper.run_one_path();
        for(idx_t i = 0; i < num_ci; i++)
            path_calib_swaption_payoff_idx[i] = lmm_calibration_helper.path_calib_swaption_payoff[i].markAsOutput(); // should be varIndex();
        

        aad_funcs->stopRecording();     
#ifdef AADC_INTERNAL_DEBUG
        aad_funcs->setVersion(aadc::AADCCompiledFunctions<mmType>::FuncVersion::BIN); // BIN);
#endif    
        auto compile_stop = std::chrono::high_resolution_clock::now();
        std::chrono::microseconds compile_time = std::chrono::duration_cast<std::chrono::microseconds>(compile_stop - compile_start);

        aad_funcs->outStats(std::cout, "LMM");

        nSIMD = sizeof(mmType) / sizeof(double);
        nRepeats=nPaths/nSIMD;
        ci_pv.resize(num_ci,0.0);

        // benchmark double case
        std::chrono::microseconds duration;
        for (int try_i = 0; try_i < 1; ++try_i)
        {
            calibration_helper<double, TestCurve<double>, LmmVol<double>> lmm_calibration_helper(data, disc_org);

            lmm_calibration_helper.xi.resize(nSteps*nFactors);

            std::vector<double> calib_prices(num_ci, 0.);

            auto start = std::chrono::high_resolution_clock::now(); 

            for (int mc_i = 0; mc_i < nPaths; ++mc_i) {
                std::copy(
                    randomSample + mc_i * lmm_calibration_helper.xi.size()
                    , randomSample + (mc_i + 1) * lmm_calibration_helper.xi.size()
                    , lmm_calibration_helper.xi.begin()
                );
                lmm_calibration_helper.run_one_path();
                for (int ci_i = 0; ci_i < num_ci; ++ci_i) {
                    calib_prices[ci_i] += lmm_calibration_helper.path_calib_swaption_payoff[ci_i];
                }
            }            
            auto stop = std::chrono::high_resolution_clock::now(); 
            duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 
            std::cout << "running MC double simulation time : " << duration.count() << std::endl;
        }
    }

    // forward only, update ci_pv
    void simulate1(const std::vector<double>& p)
    {
        std::vector<std::unique_ptr<std::thread>> threads;

        auto func = [&p, this]( const int pn_start, const int pn_end,std::vector<double>& ci_pv_chunk,
        uint64_t& ncalls_forward_th) 
        {   
            double* randomSample = (double*)(&(mmRandomSample[0]));
            std::shared_ptr<aadc::AADCWorkSpace<mmType> > ws_th(aad_funcs->createWorkSpace());    
            for(int i = 0 ; i < p.size(); i++)
                ws_th->setVal(p_idx[i],mmSetConst<mmType>(p[i]));

            
            for(idx_t pn = pn_start; pn < pn_end; pn++)
            {
                idx_t offset = pn * nSIMD * xi_idx.size();
                mmType * randoms((mmType*)(&randomSample[offset]));
                for(int i = 0; i < xi_idx.size(); i++)
                {
                    ws_th->val(xi_idx[i])=*randoms;
                    ++randoms;
                }
//                std::cout << "run fwd" << std::endl;
                aad_funcs->forward(*ws_th);
                ncalls_forward_th++;
//                std::cout << "done run fwd" << std::endl;
                for(int i = 0; i < ci_pv_chunk.size(); i++)
                {
                    for(int k = 0; k < nSIMD; k++)
                        ci_pv_chunk[i] += ws_th->val(path_calib_swaption_payoff_idx[i])[k]/double(nPaths);
                }

            }
        };

        int k(0);
        int thread_chunk(nRepeats / nThreads);
        auto start = std::chrono::high_resolution_clock::now(); 
        std::cout << std::setprecision(18);

        std::fill(ci_pv.begin(),ci_pv.end(),0.0);
        std::vector<std::vector<double>> thread_ci_pv(nThreads,ci_pv);
        std::vector<uint64_t> ncalls_forward_th(nThreads,0.0);
        for(int i=0; i< nThreads; i++) {
            threads.push_back(
                std::unique_ptr<std::thread>(
                    new std::thread(
                        func, 
                        i * thread_chunk, 
                        i == (nThreads - 1) ? nRepeats : ((i+1) * thread_chunk),
                        std::ref(thread_ci_pv[i]),
                        std::ref(ncalls_forward_th[i])
                    )
                )
            );
        }
        for(auto&& t: threads) t->join();
        auto stop = std::chrono::high_resolution_clock::now(); 
        simulate1_time = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 

        std::fill(ci_pv.begin(),ci_pv.end(),0.0);

        for(int i = 0; i < ci_pv.size(); i++)
        {
            for(int k = 0; k < nThreads; k++) {
                ci_pv[i] += thread_ci_pv[k][i];
            }
        }
        for(int k = 0; k < nThreads; k++) {
                ncalls_forward += ncalls_forward_th[k];
        }
    }

    void set_target(const std::vector<double>& p)
    {
        simulate1(p);
        target_pv = ci_pv;
        // target_pv.resize(ci_pv.size());
        // for(int i = 0; i < ci_pv.size(); i++)
        //     target_pv[i]= ci_pv[i]*2;
    }

    double cost(const std::vector<double>& p)
    {
        simulate1(p);
        double ret = 0.0;
        for(int i = 0; i < ci_pv.size(); i++)
        {
            double diff = ci_pv[i] - target_pv[i];
            ret += diff * diff;
        }
        return ret;
    }

    void simulate2(const std::vector<double>& p, bool sim1_done=false)
    {
        if (!sim1_done) simulate1(p);
        // compute cost here ...
        std::vector<double> misprice(ci_pv.size());
        for(int i = 0; i < ci_pv.size(); i++)
            misprice[i] = ci_pv[i] - target_pv[i];
        
        grad.resize(p.size());
        std::fill(grad.begin(),grad.end(),0.0);

        auto simulate_grad_thread = [&, this](
            const int pn_start, const int pn_end,
            std::vector<double>& grad_th,
            uint64_t& ncalls_forward_th,
            uint64_t& ncalls_reverse_th
        ) {
            double* randomSample = (double*)(&(mmRandomSample[0]));
            std::shared_ptr<aadc::AADCWorkSpace<mmType> > ws_th(aad_funcs->createWorkSpace());    
            for(int i = 0 ; i < p.size(); i++)
                ws_th->setVal(p_idx[i],mmSetConst<mmType>(p[i]));

            for(idx_t pn = pn_start; pn < pn_end; pn++) {
                idx_t offset = pn * nSIMD * xi_idx.size();
                mmType * randoms((mmType*)(&randomSample[offset]));
                for(int i = 0; i < xi_idx.size(); i++)
                {
                    ws_th->val(xi_idx[i])=*randoms;
                    ++randoms;
                }

                //std::cout << "run fwd" << std::endl;
                (*aad_funcs).forward(*ws_th);
                ncalls_forward_th++;
                //std::cout << "done run fwd" << std::endl;
                ws_th->resetDiff();
                for(idx_t i = 0; i < num_ci; i++)
                    ws_th->diff(path_calib_swaption_payoff_idx[i]) = mmSetConst<mmType>(2.0*misprice[i]);

                (*aad_funcs).reverse(*ws_th);
                ncalls_reverse_th++;

                for(int k = 0; k < nSIMD; k++)
                {
                    for(int i = 0; i < p.size(); i++)
                        grad_th[i] += ws_th->diff(p_idx[i])[k];
                }
            }
        };
        int thread_chunk(nRepeats / nThreads);
        std::cout << std::setprecision(18);

        std::vector<std::vector<double>> thread_grad(nThreads,grad);
        std::vector<std::unique_ptr<std::thread>> threads;

        std::vector<uint64_t> ncalls_forward_th(nThreads,0.0);
        std::vector<uint64_t> ncalls_reverse_th(nThreads,0.0);

        for(int i=0; i< nThreads; i++) {
            threads.push_back(
                std::unique_ptr<std::thread>(
                    new std::thread(
                        simulate_grad_thread, 
                        i * thread_chunk, 
                        i == (nThreads - 1) ? nRepeats : ((i+1) * thread_chunk),
                        std::ref(thread_grad[i]),
                        std::ref(ncalls_forward_th[i]),
                        std::ref(ncalls_reverse_th[i])
                    )
                )
            );
        }
        for(auto&& t: threads) t->join();

        for(int i = 0; i < grad.size(); i++)
        {
            for(int k = 0; k < nThreads; k++) {
                grad[i] += thread_grad[k][i];
            }
        }

        for(int i = 0; i < grad.size(); i++)
            grad[i] /= double(nPaths);

        for(int k = 0; k < nThreads; k++) {
                ncalls_forward += ncalls_forward_th[k];
                ncalls_reverse += ncalls_reverse_th[k];
        }

    }

    double cost(const std::vector<double>& p, std::vector<double>& grad_)
    {
        simulate2(p);
        grad_.resize(p.size());
        std::copy(grad.begin(),grad.end(),grad_.begin());
        double ret = 0.0;
        for(int i = 0; i < ci_pv.size(); i++)
        {
            double diff = ci_pv[i] - target_pv[i];
            ret += diff * diff;
        }
        return ret;
    }
    std::vector<double> params;

    // for calibrator: map unbounded parameters x to meaningful intervals
    std::vector<double> params_from;
    std::vector<double> params_to;

    double F      (double x, double a, double b) { return a + (b-a)/(1.0+std::exp(-x));}
    double F_prime(double x, double a, double b) { 
        double e = std::exp(-x);
        double d = (1+e);
        return (b-a)*e/(d*d);
    }
    double F_inv  (double p, double a, double b) { return std::log((p-a)/(b-p));}

    double operator()(const VectorXd& x, VectorXd& gradient)
    {
        // auto start = std::chrono::high_resolution_clock::now(); 
        params.resize(x.size());
        
        for(int i = 0; i < params.size(); i++)
            params[i] = F(x[i],params_from[i],params_to[i]);

        simulate2(params);
        
        // large factor to make initial step safely small
        
        for(int i = 0; i < params.size(); i++)
            gradient[i] = grad[i] * factor * F_prime(x[i],params_from[i],params_to[i]);

        double ret = 0.0;

        for(int i = 0; i < ci_pv.size(); i++)
        {
            double misprice = ci_pv[i] - target_pv[i];
            ret += misprice * misprice; // discrepancy
        }
        // auto stop = std::chrono::high_resolution_clock::now(); 
        // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 

        // std::cout << "Total squared errors " << ret << " params ";
        // for( int i = 0; i < params.size(); i++)
        //     std::cout << params[i] << " ";
        // std::cout    << " time this iteration : " << duration.count() << " time fwd+rev " << (duration.count() - simulate1_time.count()) << std::endl; 
        //std::cout.flush();
        return ret * factor;
    }
    /// operator that computes fx in all points given
    double operator()( VectorXd &x, VectorXd& grad1, double &step, const VectorXd &drt, const VectorXd &xp)
    {   
        if (!onlyF) return operator()(x,grad1);

        for(int i = 0; i < params.size(); i++)
            params[i] = F(x[i],params_from[i],params_to[i]);

        simulate1(params);

        double fx=0.0;
        for(int i = 0; i < ci_pv.size(); i++){
            double misprice = ci_pv[i] - target_pv[i];
            fx += misprice * misprice; // discrepancy
        }
        
        return fx*factor;
    }
    /// function that calculates the 1D-derivative because most iterations there is no need for the full gradient
    double getDg(const VectorXd &drt, const VectorXd &x , const double &fx0, const VectorXd grad1)  
    {   
        if (!onlyF) return grad1.dot(drt);

        auto ci_pv1 = ci_pv;

        VectorXd x2;
        x2.noalias() = x + delta * drt;  

        double fx1 = fx0/factor;
    
        std::vector<double> params2 = params;
        double fx2 = 0;

        for(int i = 0; i < params.size(); i++)
            params2[i] = F(x2[i],params_from[i],params_to[i]);

        simulate1(params2);

        for(int i = 0; i < ci_pv.size(); i++){
            double misprice = ci_pv[i] - target_pv[i];
            fx2 += misprice * misprice; // discrepancy
        }

        double res = (fx2 - fx1) /delta * factor;
        ci_pv=ci_pv1;
        return res; 
    }

    // function to get the full gradient when needed
    void getGrad(const VectorXd& x, VectorXd& gradient)
    {   
        if (!onlyF) return;
        params.resize(x.size());
        
        for(int i = 0; i < params.size(); i++)
            params[i] = F(x[i],params_from[i],params_to[i]);

        // true because the simulate1 results are already calculated in the normal operator()
        simulate2(params,true); 
        
        // large factor to make initial step safely small
        for(int i = 0; i < params.size(); i++)
            gradient[i] = grad[i] * factor * F_prime(x[i],params_from[i],params_to[i]);

    }
    const int getNcallsR(){
        return ncalls_reverse;
    }
    const int getNcallsF(){
        return ncalls_forward;
    }
    const void setDelta(double d) {
        delta=d;
    }
    const double getDelta() {
        return delta;
    }
    const void setOnlyF(bool b) {
        onlyF=b;
    }
    const bool getOnlyF() {
        return onlyF;
    }
    
};

template<class mmType>
int run_example(int num_threads, bool onlyF)
{
    auto data = test_data();
    TestCurve<idouble>  disc;
    TestCurve<double>  disc_org;

    calibration_evaluator<mmType> calibrator(data, disc, disc_org, num_threads);
    calibrator.setOnlyF(onlyF);
    LmmVol<double> vol;
    vol.setStartingPoint0();

    //if (argc == 1) return 1;

    calibrator.params_from.resize(vol.p.size());
    calibrator.params_to.resize(vol.p.size());

    calibrator.params_from = vol.params_from();
    calibrator.params_to   = vol.params_to();

    VectorXd x(vol.p.size());
    for(int i = 0; i < vol.p.size(); i++)
    {
        // x[i] = calibrator.F_inv(0.75*vol.p[i],calibrator.params_from[i],calibrator.params_to[i]);
        x[i] = calibrator.F_inv(0.9*vol.p[i],calibrator.params_from[i],calibrator.params_to[i]);
    }
    // std::cout << "x = \n" << x.transpose() << std::endl;


    calibrator.set_target(vol.p);

    LBFGSParam<double> param;
    param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
    // param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
    param.epsilon = 1e-6;
    // param.max_iterations = 100; 
    param.max_linesearch = 200;

    // Create solver and function object
    LBFGSSolver<double,LineSearchBacktracking> solver(param);
    // LBFGSSolver<double,LineSearchBracketing> solver(param);
    

    double cost;
    
    auto start = std::chrono::high_resolution_clock::now(); 

    auto niters = solver.minimize(calibrator, x, cost);

    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << std::endl;
    std::cout << std::get<0>(niters) << " iterations, ";
    std::cout << std::get<1>(niters) << " iterations LS, ";
    std::cout << "ncallsR: " << calibrator.getNcallsR() << ", ncallsF: " << calibrator.getNcallsF() << std::endl;
    std::cout << duration.count() << " miliseconds" << std::endl;
    // std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "cost(x) = " << cost << std::endl;

    return 0;

}

int main(int argc, char* argv[])
{
    // int num_threads = 1;
    int avx_type = 512;
    int num_threads = 20;

    if (argc > 1) {
        avx_type = atoi(argv[1]) == 512 ? 512 : 256;
    }
    if (argc > 2) {
        num_threads = atoi(argv[2]);
    }

#if AADC_512 
    if (avx_type == 512) {
        std::cout << "AVX512" << " Num Threads : " << num_threads << " Before:" << std::endl;
        run_example<__m512d>(num_threads,false);
        std::cout << "\nAVX512" << " Num Threads : " << num_threads << " After:" << std::endl;
        return run_example<__m512d>(num_threads,true);
    }
#endif
    std::cout << "AVX2" << " Num Threads : " << num_threads << " Before:" << std::endl;
    run_example<__m256d>(num_threads,false);
    std::cout << "\nAVX2" << " Num Threads : " << num_threads << " After:" << std::endl;
    return run_example<__m256d>(num_threads,true);
}
