// #include "stdafx.h"

#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <ctime>

#include <iomanip>
#include <sstream>
#include <math.h>

#include <immintrin.h>

#include <unordered_map>

#include <memory>
#include <cstdlib>

#include "runtime.h"
#include <aadc/idouble.h>

class tdouble {
public:
    tdouble(const double& a) : val(a) {
        std::cout << "tdouble(double)" << std::endl;
    }
    tdouble(const tdouble& o) : val(o.val) {
        std::cout << "tdouble(tdouble)" << std::endl;
    }
    tdouble() {}
    ~tdouble() {
        std::cout << "~tdouble()" << std::endl;
    }
    tdouble& operator = (const tdouble& o) {
        std::cout << "operator = (tdouble)" << std::endl;
        val = o.val;
        return *this;
    }

    double val;
};


int main()
{

    std::vector<tdouble> a;

    for (int i  = 0; i < 1000; ++i) {
        std::cout << "Elem : " << i;
        a.push_back(tdouble(i*0.1));
    }

    return 0;
}

