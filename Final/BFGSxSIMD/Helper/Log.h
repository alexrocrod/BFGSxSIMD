/// Alexandre Rodrigues , 29/1/2021 

#pragma once

#include <iostream>
#include <Eigen/Core>
#include <chrono>

namespace Log{
    /// log function to print the result of x and fx 
    void logX(const Eigen::VectorXd& x, const double& fx, bool log_x){
        if (!log_x) return;
        std::cout << std::endl;
        std::cout << "x = \n" << x.transpose() << std::endl;
        std::cout << "f(x) = " << fx << std::endl;
    }

    /// Timer class for simple timing of any code, returning in milli or microseconds
    class Timer{
    private:
        std::chrono::time_point<std::chrono::high_resolution_clock> startpoint;
    public:
        Timer() {
            startpoint = std::chrono::high_resolution_clock::now();
        }
        double stopMicro(){
            auto endpoint =std::chrono::high_resolution_clock::now();
            auto start = std::chrono::time_point_cast<std::chrono::microseconds>(startpoint).time_since_epoch().count();
            auto end = std::chrono::time_point_cast<std::chrono::microseconds>(endpoint).time_since_epoch().count();
            auto duration = end - start;
            return duration ;
        }
        double stop(){
            auto endpoint =std::chrono::high_resolution_clock::now();
            auto start = std::chrono::time_point_cast<std::chrono::milliseconds>(startpoint).time_since_epoch().count();
            auto end = std::chrono::time_point_cast<std::chrono::milliseconds>(endpoint).time_since_epoch().count();
            auto duration = end - start;
            return duration ;
        }
    };
}
