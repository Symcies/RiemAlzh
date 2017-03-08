### PROJECT

cmake_minimum_required(VERSION 3.2)
project(Longitudina CXX)
Set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

option(BUILD_TEST "BUILD TEST FOR RiemAlzh" ON)
option(BUILD_BENCH "BUILD BENCHMARKS FOR RiemAlzh" OFF)

if (APPLE)
    set(CMAKE_MACOSX_RPATH 1)
endif()

### ARMADILLO

find_Package(Armadillo 4.32 QUIET)
if (ARMADILLO_FOUND)
    Include_Directories(${ARMADILLO_INCLUDE_DIRS})
    Message("Armadillo has been found")
else()
    Message(FATAL_ERROR " -- IMPOSSIBLE TO FIND ARMADILLO -- ")
endif()


### Set new project folders / files

Set(SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}

        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/algorithm
        ${CMAKE_CURRENT_SOURCE_DIR}/src/manifolds
        ${CMAKE_CURRENT_SOURCE_DIR}/src/models
        ${CMAKE_CURRENT_SOURCE_DIR}/src/models/network
        ${CMAKE_CURRENT_SOURCE_DIR}/src/models/scalar
        ${CMAKE_CURRENT_SOURCE_DIR}/src/manifolds/base_manifold
        ${CMAKE_CURRENT_SOURCE_DIR}/src/parameters
        ${CMAKE_CURRENT_SOURCE_DIR}/src/random_variables
        ${CMAKE_CURRENT_SOURCE_DIR}/src/samplers
        ${CMAKE_CURRENT_SOURCE_DIR}/src/io
        ${CMAKE_CURRENT_SOURCE_DIR}/src/realizations
        ${CMAKE_CURRENT_SOURCE_DIR}/src/io/outputs
        ${CMAKE_CURRENT_SOURCE_DIR}/src/observations
        ${CMAKE_CURRENT_SOURCE_DIR}/src/support
        ${CMAKE_CURRENT_SOURCE_DIR}/src/support/linear_algebra
        ${CMAKE_CURRENT_SOURCE_DIR}/src/support/linear_algebra/armadillo

#TODO: Compile test apart from main code.
#Are the current test files utils instead of tests?
        ${CMAKE_CURRENT_SOURCE_DIR}/test
        )

#TODO: Change the recurse and include all files manually
file( GLOB_RECURSE SRC src/* )
file( GLOB_RECURSE TEST test/* )
file( GLOB_RECURSE EXAMPLES examples/* )

### TINYXML2

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lib/tinyxml2)

### LIBRARIES

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

include_directories(${SOURCE_FILES})

# Creation of the src library, linked with armadillo and tinyxml
Add_Library(riem-lib
  STATIC
  ${SRC}
#TODO: Compile test apart from main code.
  ${TEST})
Target_Link_Libraries(riem-lib
  ${ARMADILLO_LIBRARIES}
  tinyxml2)

### EXECUTABLE

#Creation of the executable, linked with the src library
add_executable(Longitudina
        ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)

Target_Link_Libraries(Longitudina riem-lib)

### TESTS

if(BUILD_TESTS)
    enable_tests()
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lib/googletest)
    #add_subdirectory(test)
endif()
