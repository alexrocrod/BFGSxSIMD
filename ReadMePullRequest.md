This readme explains the initial modifications done by Alexandre Rodrigues <<xandre.rodrigues@hotmail.com>>.

These alterations divide the interface interaction in 3 different functions. 

We based this approach on the fact that the full gradient is only needed in the very end of the line search, during the various iterations of line search we only need the 1D derivative (except for Armijo condition that never needs this calculation).

The normal `operator()` is substituted by a new operator that calculates only the function at that point.

The function getDg() calculates the derivative by finite difference: calculates the function at the point x+delta and then: f'(x) = [ f(x + &Delta; ) - f(x) ] / &Delta;. We used `delta` as `1e-7` but you should check if there is another value that gives you a better result, in some cases of our tests we got a worse result but much faster runtime with `delta = 1e-5`.

After exiting the line search we need to calculate the full gradient because it's needed for the correction applied in the main LBFGS algorithm. This is done by calling the function `getGrad` of the interface, that only computes the gradient (normally part of the `operator()`). In some cases it’s possible to get more performance if the results of the new `operator()` can be used in the gradient calculation.

## An Example 

For the [Rosenbrock function](https://en.wikipedia.org/wiki/Rosenbrock_function):

```cpp
#pragma once

#include <Eigen/Core>
#include <LBFGS.h>
#include <vector>

using Eigen::VectorXd;

class RosenbrockNew
{
private:
    int ncalls = 0;
    int n;
    const double delta = 1e-7;

public:
    RosenbrockNew (int ni) : n(ni) {}

    const int getNcalls() {
        return ncalls;
    }
    // traditional double operator
    double operatorDouble(const VectorXd& x, VectorXd& grad)
    {   
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
    /// new operator() 
    double operator()( const VectorXd &x, VectorXd& grad )
    {   
        ncalls++;

        if (ncalls==1) return operatorDouble(x,grad);

        return getForward(x);
    }
    /// function that calculates the 1D-derivative because most iterations there is no need for the full gradient
    double getDg( const VectorXd &drt, const VectorXd &xp , const double &fx0)   
    {        
        VectorXd x2;
        x2.noalias() = xp + delta * drt;    
    
        double fx2 = getForward(x2);

        double res = (fx2 - fx0) / delta;

        return res; 
    }
    /// function that only calculates the function value
    double getForward( const VectorXd &x)
    {   
        double fx = 0.0;
        for(int i = 0; i < n ; i += 2){
            double t1 = 1.0 - x[i];
            double t2 = 10 * (x[i + 1] - x[i] * x[i]);
            fx += t1 * t1 + t2 * t2;
        }  
        return fx;   
    }
    // function to get the full gradient when needed
    void getGrad(const VectorXd& x, VectorXd& grad)
    {   
        for(int i = 0; i < n ; i += 2){
            double t1 = 1.0 - x[i];
            double t2 = 10 * (x[i + 1] - x[i] * x[i]);
            grad[i + 1] = 20 * t2;
            grad[i]     = -2.0 * (x[i] * grad[i + 1] + t1);
        }  
    }
};
```
1\.	New `operator()`: same usage as the normal one, in the first calculation we need the gradient but afterwards only need the objective function value. So works the same way but the gradient is only updated on the 1<sup>st</sup> call.

```cpp
/// new operator() 
double operator()( const VectorXd &x, VectorXd& grad )
{   
    ncalls++;

    if (ncalls==1) return operatorDouble(x,grad);

    return getForward(x);
}
double operatorDouble(const VectorXd& x, VectorXd& grad)
{   
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
double getForward( const VectorXd &x)
{   
    double fx = 0.0;
    for(int i = 0; i < n ; i += 2){
        double t1 = 1.0 - x[i];
        double t2 = 10 * (x[i + 1] - x[i] * x[i]);
        fx += t1 * t1 + t2 * t2;
    }  
    return fx;   
}
```

2\.	getDg():

 -	Input: direction (`drt`), point (`x`), function value on the previous call of the `operator()`(`fx1`);

 -	Output: returns the result of the 1D derivative;

 -	Algorithm: calculation of the derivative with finite differences, example:
    
    -	`Delta = 1e-7`, should be adjusted for the application;

    -   `x1= x` ;
    -   `x2 = x + delta * drt` ;

    -	Calculation of the function value in `x2`:`fx2` ;

    -	Return `(fx2-fx1) / delta` ;


```cpp
    /// function that calculates the 1D-derivative because most iterations there is no need for the full gradient
    double getDg( const VectorXd &drt, const VectorXd &xp , const double &fx0)   
    {        
        VectorXd x2;
        x2.noalias() = xp + delta * drt;    
    
        double fx2 = getForward(x2);

        double res = (fx2 - fx0) / delta;

        return res; 
    }
```
3)	getGrad():

-	Input: point (`x`);
-	Output: gradient passed by reference (`grad`);
-	Algorithm: calculation of the gradient, in some cases we can get better performance by using the results of the `operator()` calculation and saving 1 call of the forward algorithm.

```cpp
    // function to get the full gradient when needed
    void getGrad(const VectorXd& x, VectorXd& grad)
    {   
        for(int i = 0; i < n ; i += 2){
            double t1 = 1.0 - x[i];
            double t2 = 10 * (x[i + 1] - x[i] * x[i]);
            grad[i + 1] = 20 * t2;
            grad[i]     = -2.0 * (x[i] * grad[i + 1] + t1);
        }  
    }
```
## Changes in Line Search algorithm: 

Changes in LineSearchBacktraking.h, only 1 line changed (`getDg`):
```cpp
(...) at line 82:
                //const Scalar dg = grad.dot(drt);
                const Scalar dg = f.getDg(drt, xp, fx);                
(...)
```

Changes needed for Bracketing and Nocedal Wright algorithms are similar, box-constrained algorithm was not used.

## Changes in LBFGS algorithm: 

Changes in LBFGS.h (only after the beginning of the for loop), only 1 line added (`getDg`):
```cpp
(...)
    // Number of iterations used
    int k = 1;
    for( ; ; )
    {
        // Save the curent x and gradient
        m_xp.noalias() = x;
        m_gradp.noalias() = m_grad;

        // Line search to update x, fx and gradient
        LineSearch<Scalar>::LineSearch(f, fx, x, m_grad, step, m_drt, m_xp, m_param);
        
        // needs the full gradient
        f.getGrad(x,m_grad); // << added

        // New gradient norm
        gnorm = m_grad.norm();
(...)
```

## Conclusion

These changes give us a new way of interaction with the BFGS method that can, in some specific cases, give approximately equal results with less computational cost. 

The best results are on applications that have very complex function value and gradient calculations. One should also try to take advantage of the function value calculated in the `operator()` in the calls of `getDg` and `getGrad`, this is the best and easier way to save computational resources.

### Results

With only these changes, the number of iterations doesn't change, and the processing time is worse in this specific example, takes 5% to 15% more time.

Using the Strong Wolfe condition, we got bad results with Rosenbrock, the number of calls increased very much (by 20% in Nocedal Wright and by 50x in the others). One should test for the specific application. These results were only achieved by using delta as `1e-2`, smaller or larger delta throws the error of max line search iterations.

We also worked in a much more computational intensive application, using the function value calculated in the `operator()` we saved some time in  `getDg` and `getGrad`. With the Armijo condition (1D derivative not needed) we got time savings between 15% and 20%. 

### Disadvantages

Tests done with this Rosenbrock give some problems such as `the moving direction increases the objective function value` when using random initial points.

We believe this happen only in this type of functions and because the random initial point is not a good point for the derivative calculation.
Changing the delta value can help depending on the application.

### TO-DO:

- Check if the parameters of LBFGS can be tuned to give better results specific for this implementation;

- Improvement of this work integration on the Nocedal Wright Algorithm;

- Adapt this work for LBFGS-B algorithm if possible.

This is a work in progess, feel free to give us feedback or improvements.

