#include <Rosenbrock.h>
#include <BfgsHelper.h>
#include <Log.h>
#include <iostream>
#include <aadc/aadc.h>

using Eigen::VectorXd;
using namespace LBFGSpp;
using namespace aadc;
using Log::Timer;
using Log::logX;

typedef __m512d mmType;
// typedef __m256d mmType;

int main()
{   
    // Parameters definition 
    LBFGSParam<double> param;
    param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
    // param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
    // param.max_linesearch = 200; 
    LBFGSSolver<double,LineSearchBacktracking> solver_double(param);                               

    LBFGSSolver<double,LineSearchBacktracking> solver_nopoly(param);                                                     

    LBFGSSolver<double,LineSearchBacktracking> solver_poly(param);                               
    
    const int tests_per_n = 1024;    // default: 1024
    const int max_n = 100;          // default: 24
    const float tol = 1e-4;         // default: 1e-4
    Timer timerAll;
    float ms_double(0), ms_nopoly(0), ms_poly(0);
    
    for( int n=10; n <= max_n; n += 2 ) 
    {
        std::cout << "n = " << n << std::endl;

        std::cout.setstate(std::ios_base::failbit); // override cout to disable evaluation version warning
         
        int niter[3] = {0};
        int niterLS[3] = {0};
        VectorXd x, x0;

        Rosenbrock fun(n); // old LBFGS++

        // Without polyfit
        std::shared_ptr<aadc::AADCFunctions<mmType> > aad_funcs2; 
        std::shared_ptr<aadc::AADCWorkSpace<mmType> > ws2;
        std::vector<aadc::AADCArgument> v_index2;
        std::vector<aadc::AADCResult> f_res2;
        {
            std::vector<idouble> x(n);
            idouble fx(0.0), t1, t2;

            aad_funcs2 = std::make_shared<aadc::AADCFunctions<mmType> >();

            aad_funcs2->startRecording();

                for(int i = 0; i < n ; i++ ){
                    v_index2.push_back(x[i].markAsInput());
                }
                for(int i = 0; i < n; i += 2){
                    t1 = 1.0 - x[i];
                    t2 = 10 * (x[i + 1] - x[i] * x[i]);
                    fx += t1 * t1 + t2 * t2;
                }
                f_res2.push_back(fx.markAsOutput()); 

            aad_funcs2->stopRecording();

            ws2 = aad_funcs2->createWorkSpace();
        }

        BfgsHelper<mmType> fun_nopoly(aad_funcs2, ws2, v_index2,f_res2);
        fun_nopoly.setPolyFitOrder(0);

        // With polyfit
        std::shared_ptr<aadc::AADCFunctions<mmType> > aad_funcs; 
        std::shared_ptr<aadc::AADCWorkSpace<mmType> > ws;
        std::vector<aadc::AADCArgument> v_index;
        std::vector<aadc::AADCResult> f_res;
        {
            std::vector<idouble> x(n);
            idouble fx(0.0), t1, t2;
            aad_funcs = std::make_shared<aadc::AADCFunctions<mmType> >();
            aad_funcs->startRecording();
                for(int i = 0; i < n ; i++ ){
                    v_index.push_back(x[i].markAsInput());
                }
                for(int i = 0; i < n; i += 2){
                    t1 = 1.0 - x[i];
                    t2 = 10 * (x[i + 1] - x[i] * x[i]);
                    fx += t1 * t1 + t2 * t2;
                }
                f_res.push_back(fx.markAsOutput()); 
            aad_funcs->stopRecording();
            ws = aad_funcs->createWorkSpace();
        }   
        BfgsHelper<mmType> fun_poly(aad_funcs, ws, v_index, f_res);

        std::cout.clear(); // re-enable cout normal behaviour
   

        for( int test=0; test < tests_per_n; test++ ) 
            {
                x0 = VectorXd::Random(n)*0.2 + VectorXd::Constant(n,1.0); //random initial x from 0.8 to 1.2
                int niteri[3] = {0};
                int niterLSi[3] = {0};
                double fx;

                // old
                Timer timerDouble;
                x = x0; int i = 0;
                
                std::tie(niteri[i], niterLSi[i]) = solver_double.minimize(fun, x, fx);
                
                niter[i] +=  niteri[i]; niterLS[i] +=  niterLSi[i];
                assert(((x.array() - 1.0).abs() < tol ).all());

                ms_double += timerDouble.stopMicro();
                
                // without polyfit
                Timer timer_nopoly;
                x = x0; i++;
                
                std::tie(niteri[i], niterLSi[i]) = solver_nopoly.minimize(fun_nopoly, x, fx);
                
                niter[i] +=  niteri[i]; niterLS[i] +=  niterLSi[i];
                assert(((x.array() - 1.0).abs() < tol ).all());
                
                ms_nopoly += timer_nopoly.stopMicro();

                // with polyfit
                Timer timer_poly;
                x = x0; i++;
                
                std::tie(niteri[i], niterLSi[i]) = solver_poly.minimize(fun_poly, x, fx);
                
                niter[i] +=  niteri[i]; niterLS[i] +=  niterLSi[i];
                assert(((x.array() - 1.0).abs() < tol ).all());
                
                ms_poly += timer_poly.stopMicro();
            }
            std::cout << "  Old        : " << (niter[0] / tests_per_n) << " iterations, ";
            std::cout << (niterLS[0] / tests_per_n) << " iterations LS, "  << (fun.getNcalls() / tests_per_n) << " calls, " << ms_double << " micros"<< std::endl;
            std::cout << "  No Polyfit : " << (niter[1] / tests_per_n) << " iterations, ";
            std::cout << (niterLS[1] / tests_per_n) << " iterations LS, "  << (fun_nopoly.getNcallsR() / tests_per_n) << " callsR, " << (fun_nopoly.getNcallsF() / tests_per_n) << " callsF, "<< ms_nopoly << " micros"<< std::endl;
            std::cout << "  Polyfit    : " << (niter[2] / tests_per_n) << " iterations, ";
            std::cout << (niterLS[2] / tests_per_n) << " iterations LS, "  << (fun_poly.getNcallsR() / tests_per_n) << " callsR, " << (fun_poly.getNcallsF() / tests_per_n) << " callsF, "<< ms_poly << " micros"<< std::endl;
            
        }
    std::cout << std::endl << "Total duration: "<< timerAll.stop()  <<" ms" << std::endl;  
    return 0;
}