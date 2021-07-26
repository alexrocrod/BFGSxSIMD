#include <Eigen/Core>
#include <iostream>
#include <LBFGS.h>

using Eigen::VectorXf;
using Eigen::MatrixXf;
using namespace LBFGSpp;

class Rosenbrock
{
private:
    int n;
    uint64_t ncalls = 0;
    float delta = 1e-7;
public:
    Rosenbrock(int n_) : n(n_) {}
    float operator()(const VectorXf& x, VectorXf& grad)
    {
        ncalls++;
        if (ncalls != 1) return getForward(x);

        float fx = 0.0;
        for(int i = 0; i < n; i += 2)
        {
            float t1 = 1.0 - x[i];
            float t2 = 10 * (x[i + 1] - x[i] * x[i]);
            grad[i + 1] = 20 * t2;
            grad[i]     = -2.0 * (x[i] * grad[i + 1] + t1);
            fx += t1 * t1 + t2 * t2;
        }
        return fx;
    }
    float getForward(const VectorXf& x)
    {
        float fx = 0.0;
        for(int i = 0; i < n; i += 2)
        {
            float t1 = 1.0 - x[i];
            float t2 = 10 * (x[i + 1] - x[i] * x[i]);
            fx += t1 * t1 + t2 * t2;
        }
        return fx;
    }
    float getDg(const VectorXf& drt, const VectorXf& x1, const float& fx)
    {
        float res = 0;
        VectorXf x2 = x1 + delta * drt;
        float fx2 = getForward(x2);
        res = ( fx2 - fx ) / delta;
        return res;
    }
    void getGrad(const VectorXf& x, VectorXf& grad)
    {
        for(int i = 0; i < n; i += 2)
        {
            float t1 = 1.0 - x[i];
            float t2 = 10 * (x[i + 1] - x[i] * x[i]);
            grad[i + 1] = 20 * t2;
            grad[i]     = -2.0 * (x[i] * grad[i + 1] + t1);
        }
    }
};

int main()
{
    const int n = 10;
    LBFGSParam<float> param;
    LBFGSSolver<float> solver(param);
    Rosenbrock fun(n);

    VectorXf x = VectorXf::Zero(n);
    float fx;
    int niter = solver.minimize(fun, x, fx);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;

    return 0;
}
