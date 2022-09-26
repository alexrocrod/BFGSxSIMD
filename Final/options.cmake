set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set (CMAKE_CXX_STANDARD 14)

option(AADC_512  "Enable support for AVX512" ON)

message("CMAKE_CURRENT_LIST_DIR is ${CMAKE_CURRENT_LIST_DIR}")

include_directories(${CMAKE_CURRENT_LIST_DIR}/aadc/include)
include_directories(${CMAKE_CURRENT_LIST_DIR}/3rdparty)
include_directories(${CMAKE_CURRENT_LIST_DIR}/3rdparty/Eigen)

if (AADC_512)
	set(AADC_HOST_ARCH "skylake-avx512")
	set(AADC_AVX_PATH "avx512")
	add_definitions(-DAADC_512=1)
else()
	set(AADC_HOST_ARCH "haswell")
	set(AADC_AVX_PATH "avx2")
	add_definitions(-DAADC_512=0)
endif()

if (${UNIX})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=${AADC_HOST_ARCH}")
endif()

add_library(aadc STATIC IMPORTED)
SET_TARGET_PROPERTIES(aadc PROPERTIES
	IMPORTED_LOCATION ${CMAKE_CURRENT_LIST_DIR}/aadc/lib/libaadc-${AADC_AVX_PATH}.a)

find_package (Threads)

link_libraries(aadc ${CMAKE_THREAD_LIBS_INIT})
