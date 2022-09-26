AADC Demo Package
=================

PREREQUISITES:
===========

Windows:
Visual Studio 2017 or 2019

Linux:
C++ compiler. Supported compilers: g++, clang, intel c++
To build examples you need cmake.

Folder structure:
=============================
include/                - core AADC library header files.
lib/                    - core AADC library lib files.
3rdparty/               - tools for tests and benchmarks. NOT used by core AADC library
aadc-getting-started/   - examples from the manual. START HERE.
python/                 - AADC Python MODULE
prettyPrinters/         - pritty printers to aim debugging
docker/                 - dockerfile with system config examples

Suggested IDE configuration:
===========================
Install Visual Studio Code. 

Add plugins:
for C/C++   -  ms-vscode.cpptools
for cmake   -  ms-vscode.cmake-tools

You can select appropriate build configuration using CMake Tools(configured in cmake-variants.yaml)

TO BUILD FROM Visual Studio 2017/2019
=============================
1) Use solution file from visualstudio/20XX/AADC.sln
2) Using CMAKE. 
Right click at AADC package folder -> select "Open with Visual Studio". Visual Studio will use cmake integration to build and run examples

TO BUILD MANUALLY ON LINUX:
=============================
You can also check docker/test/dockerfile for dependencies and build commands

Release build:
$ mkdir build
$ cd build
$ cmake ..
$ make -j

Debug build:
$ mkdir build_debug
$ cd build_debug
$ cmake -DCMAKE_BUILD_TYPE=Debug ..
$ make -j

Release build with avx512 support:
$ mkdir build_avx512
$ cd build_avx512
$ cmake -DAADC_512=1 ..
$ make -j

You can check compile_commands.json to see exact compiler options used.

TO RUN/DEBUG AADC EXAMPLES:
=============================
See aadc-getting-started/README.txt for explanations.

build$ ./aadc-getting-started/manual 1
z:1.39121
avx[0] f(4,2,0.4)=11.0232
avx[1] f(3,2,0.3)=6.04965
avx[2] f(2,2,0.2)=3.32012
avx[3] f(1,2,0.1)=1.82212
avx[0] df/dx=5.51159,df/dy=-11.0232,df/dz=11.0232
avx[1] df/dx=3.02482,df/dy=-4.53724,df/dz=6.04965
avx[2] df/dx=1.66006,df/dy=-1.66006,df/dz=3.32012
avx[3] df/dx=0.911059,df/dy=-0.45553,df/dz=1.82212
Example 1 is done



DOCUMENTATION:
=============================
Make sure to check relevant pages of the manual for more information and
an explanation of the examples (see Manual.pdf on the left). If you have
any questions, contact us at info@matlogica.com. Once you grasp the basics of AADC,
feel free to request more advanced examples and benchmarks.


