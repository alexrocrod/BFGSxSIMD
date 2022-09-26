cmake_minimum_required(VERSION 3.8)

project(ORE CXX)

#include(CTest)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include("cmake/commonSettings.cmake")

add_library(aadc STATIC IMPORTED)

set (BUILD_FOR_MARCH haswell)

option(AADC_512  "Enable AVX 512 for AADC" ON)

if(AADC_512)
        add_definitions(-DAADC_512=1)
        set (BUILD_FOR_MARCH skylake-avx512)
endif()


set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=${BUILD_FOR_MARCH}")

SET_TARGET_PROPERTIES(aadc PROPERTIES
        IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/../../Final/aadc/lib/libaadc-avx512.a)


include_directories(../../Final/aadc/include)
include_directories(../../Final/aadc/3rdparty)

option(AADC_QL  "Enable AAD for QL" ON)
option(AADC_PASSIVE  "Enable passive AAD" OFF)

if(AADC_QL)
        add_definitions(-DAADC_QL=1)
        add_definitions(-DAADC_ALLOW_TO_PASSIVE_BOOL=1)
endif()

if(AADC_PASSIVE)
        add_definitions(-DAADC_PASSIVE)
endif()

include_directories(ThirdPartyLibs)


add_subdirectory("QuantLib")
add_subdirectory("QuantExt")
add_subdirectory("OREData")
add_subdirectory("OREAnalytics")
add_subdirectory("App")
