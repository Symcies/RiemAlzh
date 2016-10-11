cmake_minimum_required(VERSION 3.5)
project(Longitudina)



####################################
### Set new project folders / files
####################################
set(PARAMETERS
        Parameters
        Parameters/CandidateRandomVariables)

set(OUTPUTS
        Outputs
        Outputs/RandomVariableRealizations)

set(UTILITIES
        Utilities
        Utilities/MatrixFunctions)

set(SAMPLERS
        Samplers
        Samplers/AbstractSampler
        Samplers/BlockedGibbsSampler
        Samplers/HMWithinGibbsSampler
        )

set(RANDOMVARIABLES
        RandomVariables
        RandomVariables/AbstractRandomVariable
        RandomVariables/GaussianRandomVariable
        RandomVariables/LaplaceRandomVariable
        RandomVariables/ConstantRandomVariable
        )

set(MANIFOLDS
        Manifolds
        Manifolds/BaseManifold
        Manifolds/BaseManifold/AbstractBaseManifold
        Manifolds/BaseManifold/LogisticBaseManifold
        Manifolds/AbstractManifold
        Manifolds/PropagationManifold)

set(MODELS
        Models
        Models/AbstractModel
        Models/UnivariateModel
        Models/LongitudinalModel)

set(ALGORITHM
        Algorithm
        Algorithm/Algorithm)

set(TESTS 
        Tests
        Tests/TestAssert)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

set(SOURCE_FILES main.cpp)


####################################
### Set ITK
####################################
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

set(ITK_DIR /aramis/dartagnan2/Software/ITK/ITKb-4.5.1)

#find_Package(ITK REQUIRED)

if(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
else(ITK_FOUND)
    MESSAGE(FATAL ERROR "ITK not found. Please set ITK_DIR.")
endif(ITK_FOUND)



####################################
### Set GNUPLOT
####################################




####################################
### Libraries
####################################



####################################
### Add Executable
####################################
add_executable(Longitudina2
        ${SOURCE_FILES}
        ${UTILITIES}
        ${SAMPLERS}
        ${RANDOMVARIABLES}
        ${MODELS}
        ${ALGORITHM}
        ${MANIFOLDS}
        ${PARAMETERS}
        ${OUTPUTS}
        ${TESTS})
