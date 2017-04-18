include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

SET(CMAKE_BUILD_TYPE DEBUG)

file(GLOB test_files
        RunTests.cxx
        test_global.h
        test_global.cxx
        integration/*.cxx
        integration/*.h
        unit/algorithm/*.cxx
        unit/algorithm/*.h
        unit/io/*.cxx
        unit/io/*.h
        unit/models/*.cxx
        unit/models/*.h)
message("test f " ${test_files})
add_executable(RunUnitTests ${test_files})
target_link_libraries(RunUnitTests gtest gtest_main riem-lib)

add_test(run-all RunUnitTests)
