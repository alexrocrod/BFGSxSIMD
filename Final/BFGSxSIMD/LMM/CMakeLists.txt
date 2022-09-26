cmake_minimum_required(VERSION 3.10)

add_executable(lmm_test lmm_test.cpp lmm_test_data.cpp)

target_link_libraries (lmm_test ${CMAKE_THREAD_LIBS_INIT})
