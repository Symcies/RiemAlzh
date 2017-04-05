### PROJECT

cmake_minimum_required(VERSION 3.2)
project(Longitudina CXX)

option(BUILD_TESTS "Build the tests ?" ON)


##This is a C++11 project!
include(CheckCXXCompilerFlag)
Set(ENABLE_CXXFLAGS_TO_CHECK
        -std=gnu++1z
        -std=c++1z
        -std=gnu++1y
        -std=c++1y
        -std=gnu++11
        -std=c++11
        -std=gnu++0x
        -std=c++0x)
Set(CXX_COMPILER_ENABLED FALSE)
foreach(flag ${ENABLE_CXXFLAGS_TO_CHECK})
    string(REPLACE "-std=" "_" flag_var ${flag})
    string(REPLACE "+" "x" flag_var ${flag_var})
    CHECK_CXX_COMPILER_FLAG("${flag}" COMPILER_HAS_CXX_FLAG${flag_var})
    if(COMPILER_HAS_CXX_FLAG${flag_var})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
        Set(CXX_COMPILER_ENABLED TRUE)
        break()
    endif()
endforeach()
if (NOT CXX_COMPILER_ENABLED)
    PRINT_ENV(${CMAKE_DEBUG})
    Message(FATAL_ERROR "The current CPP compiler does not support C++11: Please provide a good one.")
endif()

if (APPLE)
    set(CMAKE_MACOSX_RPATH 1)
endif()



### ARMADILLO
find_Package(Armadillo)
if (NOT ARMADILLO_FOUND)
    PRINT_ENV(${CMAKE_DEBUG})
    Message(FATAL_ERROR "ARMADILLO library is not found")
endif()
if (NOT ARMADILLO_VERSION_MAJOR OR "${ARMADILLO_VERSION_MAJOR}" LESS 6)
    PRINT_ENV(${CMAKE_DEBUG})
    Message(FATAL_ERROR "ARMADILLO version must be great or equal to version 6.0")
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
        ${CMAKE_CURRENT_SOURCE_DIR}/src/io/data_settings
        ${CMAKE_CURRENT_SOURCE_DIR}/src/observations
        ${CMAKE_CURRENT_SOURCE_DIR}/src/support
        ${CMAKE_CURRENT_SOURCE_DIR}/src/support/linear_algebra
        ${CMAKE_CURRENT_SOURCE_DIR}/src/support/linear_algebra/armadillo
        )

#TODO: Change the recurse and include all files manually
file( GLOB_RECURSE SRC src/* )
file( GLOB_RECURSE TEST test/* )
file( GLOB_RECURSE EXAMPLES examples/* )

####################################################################################################
#Next lines are just to work with one moddel
file(GLOB to_remove
        ${CMAKE_CURRENT_SOURCE_DIR}/src/models/network/*
        ${CMAKE_CURRENT_SOURCE_DIR}/src/models/scalar/MixedEffectsModel.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/models/scalar/MixedEffectsModel.h
        )

list(REMOVE_ITEM SRC ${to_remove})
####################################################################################################


### TINYXML2

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lib/tinyxml2)

### LIBRARIES

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  "${CMAKE_CURRENT_SOURCE_DIR}/bin")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

include_directories(${SOURCE_FILES})

# Creation of the src library, linked with armadillo and tinyxml
Add_Library(riem-lib
  STATIC
  ${SRC})
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
    enable_testing()
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lib/googletest)
    add_subdirectory(test)
endif()
