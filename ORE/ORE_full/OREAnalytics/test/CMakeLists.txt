# cpp files, this list is maintained manually

set(OREAnalytics-Test_SRC aggregationscenariodata.cpp
cube.cpp
observationmode.cpp
scenariogenerator.cpp
scenariosimmarket.cpp
sensitivityaggregator.cpp
sensitivityanalysis.cpp
sensitivityanalysisanalytic.cpp
sensitivityperformance.cpp
shiftscenariogenerator.cpp
stresstest.cpp
swapperformance.cpp
testmarket.cpp
testportfolio.cpp
testsuite.cpp)

find_package (Threads)

add_executable(orea-test-suite ${OREAnalytics-Test_SRC})
target_link_libraries(orea-test-suite ${QL_LIB_NAME})
target_link_libraries(orea-test-suite ${QLE_LIB_NAME})
target_link_libraries(orea-test-suite ${ORED_LIB_NAME})
target_link_libraries(orea-test-suite ${OREA_LIB_NAME})
target_link_libraries(orea-test-suite ${Boost_LIBRARIES})

if(AADC_QL)
    target_link_libraries(orea-test-suite aadc)
    target_link_libraries(orea-test-suite ${CMAKE_THREAD_LIBS_INIT})
endif()

add_test(NAME orea-test-suite COMMAND orea-test-suite WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
