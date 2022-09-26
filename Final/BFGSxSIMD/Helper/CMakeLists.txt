cmake_minimum_required(VERSION 3.10)

add_library(rosenbrock INTERFACE)
target_include_directories(rosenbrock INTERFACE .)
target_link_libraries (rosenbrock INTERFACE  ${CMAKE_THREAD_LIBS_INIT})


add_library(log INTERFACE)
target_include_directories(log INTERFACE .)
target_link_libraries (log INTERFACE  ${CMAKE_THREAD_LIBS_INIT})

add_executable(testPoly test.cpp)
target_link_libraries (testPoly  ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries (testPoly  rosenbrock)
target_link_libraries (testPoly  log)

add_executable(testPolyforcomp test_comparison.cpp)
target_link_libraries (testPolyforcomp  ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries (testPolyforcomp  rosenbrock)
target_link_libraries (testPolyforcomp  log)
