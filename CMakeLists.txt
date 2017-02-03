### PROJECT

cmake_minimum_required(VERSION 3.2)
project(Longitudina CXX)
Set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

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


### OPENMP

find_package(OpenMP)
if (OPENMP_FOUND)
    #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    Message("OpenMP has been found")
else()
    Message(WARNING "NOOOOO")
    MESSAGE(Warning OpenMP)    
    Set(CMAKE_C_COMPILER "/usr/local/Cellar/llvm/3.9.0/bin/clang")
    Set(CMAKE_CXX_COMPILER "/usr/local/Cellar/llvm/3.9.0/bin/clang++")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -l/usr/local/Cellar/llvm/3.9.0/lib/libomp.dylib")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /usr/local/Cellar/llvm/3.9.0/lib/libomp.dylib")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

endif()



### Set new project folders / files

Set(SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}
        
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/algorithm
        ${CMAKE_CURRENT_SOURCE_DIR}/src/manifolds
        ${CMAKE_CURRENT_SOURCE_DIR}/src/models
        ${CMAKE_CURRENT_SOURCE_DIR}/src/manifolds/base_manifold
        ${CMAKE_CURRENT_SOURCE_DIR}/src/parameters
        ${CMAKE_CURRENT_SOURCE_DIR}/src/random_variables
        ${CMAKE_CURRENT_SOURCE_DIR}/src/samplers
        ${CMAKE_CURRENT_SOURCE_DIR}/src/io
        ${CMAKE_CURRENT_SOURCE_DIR}/src/io/outputs
        
        ${CMAKE_CURRENT_SOURCE_DIR}/utilities
        ${CMAKE_CURRENT_SOURCE_DIR}/utilities/linear_algebra
        ${CMAKE_CURRENT_SOURCE_DIR}/utilities/linear_algebra/armadillo
        ${CMAKE_CURRENT_SOURCE_DIR}/utilities/matrix
        ${CMAKE_CURRENT_SOURCE_DIR}/utilities/tests
        )

file( GLOB_RECURSE SRC src/* )
file( GLOB_RECURSE UTILITIES utilities/* )

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lib/tinyxml2)

### Executable Declaration

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

add_executable(Longitudina 
        ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
        ${IO}
        ${SRC}
        ${UTILITIES})

include_directories(${SOURCE_FILES})


add_dependencies(Longitudina tinyxml2)

Target_Link_Libraries(Longitudina ${ARMADILLO_LIBRARIES})
Target_Link_Libraries(Longitudina tinyxml2)

