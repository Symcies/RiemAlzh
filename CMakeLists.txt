cmake_minimum_required(VERSION 3.5)
project(Longitudina)

set(MANIFOLDS
        Manifolds
        Manifolds/AbstractManifold
        Manifolds/MultivariateLogisticManifold)

set(MODEL
        Model
        Model/LongitudinalModel)

set(RANDOMVARIABLES
        RandomVariables
        RandomVariables/AbstractRandomVariable
        RandomVariables/GaussianRandomVariable
        RandomVariables/LaplaceRandomVariable)

set(ALGORITHM
        Algorithm
        Algorithm/Algorithm)

set(SAMPLERS
        Samplers
        Samplers/AbstractSampler
        Samplers/HastingMetropolisWithinGibbs)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

set(SOURCE_FILES
        main.cpp
        )

include_directories(
        ${MANIFOLDS}
        ${MODEL}
        ${RANDOMVARIABLES}
        ${ALGORITHM}
)

add_executable(Longitudina
        ${SOURCE_FILES}
        ${MANIFOLDS}
        ${MODEL}
        ${RANDOMVARIABLES}
        ${ALGORITHM}
        ${SAMPLERS}
        )