#export CC=/usr/bin/gcc 
#export CXX=/usr/bin/g++ 


# cpp files, this list is maintained manually
set(cpprestsdk_DIR /usr/lib/${CMAKE_LIBRARY_ARCHITECTURE}/cmake/)
find_package(cpprestsdk REQUIRED)


find_package (Threads)


include_directories(${PROJECT_SOURCE_DIR}/../aadc/3rdparty)
include_directories(../../../../Simple-WebSocket-Server)



add_executable(main main.cpp CurveCalibrate.cpp)
target_link_libraries(main ${QL_LIB_NAME})
target_link_libraries(main ${QLE_LIB_NAME})
target_link_libraries(main ${ORED_LIB_NAME})
target_link_libraries(main ${OREA_LIB_NAME})
target_link_libraries(main ${Boost_LIBRARIES})
target_link_libraries(main cpprestsdk::cpprest)


if(AADC_QL)
    target_link_libraries(main aadc)
    target_link_libraries(main ${CMAKE_THREAD_LIBS_INIT})
endif()
