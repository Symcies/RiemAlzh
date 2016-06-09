cmake_minimum_required(VERSION 3.5)
project(biomark)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++0x")

set(SOURCE_FILES main.cpp subject.cpp subject.h parameter.cpp parameter.h tools.cpp tools.h tools.cpp MCMC.cpp MCMC.h random_variable.cpp random_variable.h geodesic.cpp geodesic.h GibbsSampler.cpp GibbsSampler.h)
add_executable(biomark ${SOURCE_FILES})