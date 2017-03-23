#pragma once

typedef double ScalarType;

#include <cassert>
#include <iostream>
#include <random>
#include <unordered_map>

#include "LinearAlgebra.h"

//static std::random_device RD;
//static std::mt19937 generator(RD());

static std::mt19937 generator(1);


class AbstractRandomVariable {
public:
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
    typedef typename std::unordered_map<std::string, ScalarType > StringScalarHash;
    typedef typename std::unordered_map<int, ScalarType > IntScalarHash;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Getter(s) and Setter(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    virtual ScalarType GetParameter(std::string param_name) const = 0;
    virtual ScalarType GetParameter(int param_key) const = 0;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////
    virtual ~AbstractRandomVariable(){};
    AbstractRandomVariable(){};
    AbstractRandomVariable(const AbstractRandomVariable&){};
    virtual AbstractRandomVariable& operator=(const AbstractRandomVariable&){};

    /// Draw a new sample
    virtual double Sample() = 0;

    /// Draw multiple samples
    VectorType Samples(unsigned int samples_num);

    /// Compute the likelihood
    virtual double Likelihood(double x) = 0;

    /// Compute the log likelihood
    virtual double LogLikelihood(double x) = 0;

    /// Update the random variable parameters
    virtual void Update(StringScalarHash params) = 0;
    virtual void Update(IntScalarHash    params) = 0;

};
