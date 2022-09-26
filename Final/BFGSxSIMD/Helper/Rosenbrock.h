/// Alexandre Rodrigues , 18/2/2021 
/// traditional Rosenbrock
#pragma once

#include <Eigen/Core>
#include <LBFGS.h>
#include <vector>

using Eigen::VectorXd;

class Rosenbrock
{
private:
    uint64_t ncalls;
    int n;
    double delta = 1e-7; // default: 1e-7

public:
    Rosenbrock (int ni): 
        ncalls(0), n(ni) {}

    // traditional double operator
    double operator()(const VectorXd& x, VectorXd& grad)
    {   
        ncalls++;
        double fx = 0.0;
        for(int i = 0; i < n ; i += 2){
            double t1 = 1.0 - x[i];
            double t2 = 10 * (x[i + 1] - x[i] * x[i]);
            grad[i + 1] = 20 * t2;
            grad[i]     = -2.0 * (x[i] * grad[i + 1] + t1);
            fx += t1 * t1 + t2 * t2;
        }  
        return fx;
    }
    /// function that calculates the 1D-derivative
    double getDg(const VectorXd &drt, const VectorXd &xp , const double &fx0, VectorXd& grad)   
    {        
        return grad.dot(drt);
    }
    /// operator that computes the fx and gradient in all points given
    double operator()(VectorXd &x, VectorXd& grad, double &step, const VectorXd &drt, const VectorXd &xp  )
    {   
        return operator()(x,grad); 
    }
    // function to get the gradient when needed
    void getGrad(const VectorXd& x, VectorXd& grad)
    {   
        return; // traditional way 
    }
    const uint64_t getNcalls() {
        return ncalls;
    }
    const void setDelta(double d){
        delta = d;
    }
    const double getDelta(){
        return delta;
    }
};