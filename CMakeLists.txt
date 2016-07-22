cmake_minimum_required(VERSION 3.5)
project(Longitudina2)

set(RANDOMVARIABLES
        RandomVariables
        RandomVariables/AbstractRandomVariable
        RandomVariables/GaussianRandomVariable
        RandomVariables/LaplaceRandomVariable)

set(MANIFOLDS
        Manifolds
        Manifolds/AbstractManifold
        Manifolds/PropagationManifold)

set(MODELS
        Models
        Models/AbstractModel
        Models/LongitudinalModel)

set(ALGORITHM
        Algorithm
        Algorithm/Algorithm)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

set(SOURCE_FILES main.cpp)

add_executable(Longitudina2
        ${SOURCE_FILES}
        ${RANDOMVARIABLES}
        ${MODELS}
        ${ALGORITHM}
        ${MANIFOLDS})