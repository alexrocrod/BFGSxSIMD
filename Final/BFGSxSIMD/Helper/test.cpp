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

int main(int argc, char* argv[])
{   
    /// receives \param n as an argument (default is \c 10 )
    int n0 =10;
    int n=n0;
    if (argc > 1){
        n=atoi(argv[1]); 
        if (n == 0){
            printf("WARNING: n introduced is not greater then 1.\nUsing n as 10.\n");
            n=n0;
        }
        else if (n % 2 == 1){
            printf("WARNING: n introduced is not an even number.\nUsing n as n_intro+1.\n");
            n++;
        }
    }
    // General Function initialization
    
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

    
    BfgsHelper<mmType> fun(aad_funcs, ws, v_index, f_res);
    fun.setPolyFitOrder(0);

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

    BfgsHelper<mmType> fun_Poly(aad_funcs2, ws2, v_index2,f_res2);

    // Parameters definition 
    LBFGSParam<double> param;
    param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
    // param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
    param.max_linesearch = 20; 
    LBFGSSolver<double,LineSearchBacktracking>  solver(param);    
    LBFGSSolver<double,LineSearchBacktracking>  solver_Poly(param);  

    ///initialization
    VectorXd x0 = VectorXd::Random(n)*0.1 + VectorXd::Constant(n,1.0);  // random initial x from 0.9 to 1.1
    VectorXd x = x0;

    std::cout << "n = " << n << std::endl;
    int niter2 = 0, niter3 = 0;
    int niter_LS2 = 0 , niter_LS3 = 0;
    double fx = 0.0;

    const float tol=1e-4;   // default : 1e-4
    bool log_x=true;       // enables logging of the x and fx results

    Timer timerAll;
    
    /// specific to Rosenbrock
    
    x = x0;
    
    std::tie(niter2, niter_LS2) = solver.minimize(fun, x, fx);
    
    assert(((x.array() - 1.0).abs() < tol ).all());
    logX(x,fx,log_x);
    std::cout<<"LineSearchBack Helper : "<< niter2 << " iterations, " << niter_LS2 << " iterations LS "<< std::endl;

    /// general class

    x = x0;

    std::tie(niter3, niter_LS3) = solver_Poly.minimize(fun_Poly, x, fx);

    assert(((x.array() - 1.0).abs() < tol ).all());
    logX(x,fx,log_x);;
    std::cout<<"LineSearchBack Poly   : "<< niter3 << " iterations, " << niter_LS3 << " iterations LS "<< std::endl;

    std::cout << std::endl << "Total duration: "<< timerAll.stopMicro() <<" microsecs" << std::endl;  
    return 0;
}