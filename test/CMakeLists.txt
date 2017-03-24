include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

file(GLOB test_files RunTests.cxx usercase/*.cxx usercase/*.h)
message("test f " ${test_files})
add_executable(RunUnitTests ${test_files})
target_link_libraries(RunUnitTests gtest gtest_main riem-lib)

add_test(run-all RunUnitTests)
