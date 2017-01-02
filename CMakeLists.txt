cmake_minimum_required(VERSION 3.5)


####################################
### Project Declaration 
####################################

project(Longitudina2)



####################################
### Set new project folders / files
####################################

file( GLOB_RECURSE PARAMETERS Parameters/*)

file( GLOB_RECURSE INPUTS inputs/* )

file( GLOB_RECURSE LinearAlgebra/Utilities LinearAlgebra/Utilities/*)

file( GLOB_RECURSE LINEARALGEBRA LinearAlgebra/*)

file( GLOB_RECURSE SAMPLERS Samplers/*)

file( GLOB_RECURSE RANDOMVARIABLES RandomVariables/*)

file( GLOB_RECURSE OUTPUTS Outputs/*)

file( GLOB_RECURSE MANIFOLDS Manifolds/*)

file( GLOB_RECURSE MODELS Models/*)

file( GLOB_RECURSE ALGORITHM Algorithm/*)

file( GLOB_RECURSE TESTS Tests/* )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

set(SOURCE_FILES main.cpp)


####################################
### Set Armadillo
####################################

find_Package(Armadillo 4.32 QUIET)

if (ARMADILLO_FOUND)
    Include_Directories(${ARMADILLO_INCLUDE_DIRS})
    Message("Armadillo has been found")
else()
    Message(FATAL_ERROR " -- IMPOSSIBLE TO FIND ARMADILLO -- ")
endif()



find_package(OpenMP)
if (OPENMP_FOUND)
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    Message("OpenMP has been found")
else()
    Message("NOOOOO")    
endif()



####################################
### Set ITK
####################################

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

set(ITK_DIR /aramis/dartagnan2/Software/ITK/ITKb-4.5.1)

#find_Package(ITK REQUIRED)

#if(ITK_FOUND)
#    INCLUDE(${ITK_USE_FILE})
#    MESSAGE("ITK has been found")
#else()
#    MESSAGE(FATAL ERROR "ITK not found. Please set ITK_DIR.")
#endif(ITK_FOUND)



####################################
### Executable Declaration
####################################

add_executable(Longitudina2
        ${SOURCE_FILES}
        ${SOURCES}
        ${INPUTS}
        ${UTILITIES}
        ${SAMPLERS}
        ${RANDOMVARIABLES}
        ${MODELS}
        ${LINEARALGEBRA}
        ${ALGORITHM}
        ${MANIFOLDS}
        ${PARAMETERS}
        ${OUTPUTS}
        ${TESTS})


Target_Link_Libraries(Longitudina2
            ${ARMADILLO_LIBRARIES})
