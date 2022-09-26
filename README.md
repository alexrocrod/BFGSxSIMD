# BFGSxSIMD
Improvements and SIMD implementation for the LBFGS++ library

#### The full changes are not yet avaiable here.

We used AADC by Matlogica (matlogica.com).

Changes to LBFGS++ -> see Final/3rdparty/LBFGS.h and  Final/3rdparty/LBFGSpp/*

Helper for AADC user -> Final/BFGSxSIMD/Helper/BFGSHelper.h 

Rosenbrock test -> .cpp files in the same folder (BFGSxSIMD/Helper/), use Rosenbrock as comparison

LMM example -> BFGSxSIMD/LMM/lmm_test.cpp (changes made in calibration_evaluator class which can be disable by setting onlyF to false)

ORE example -> see ORE/howto.txt 

