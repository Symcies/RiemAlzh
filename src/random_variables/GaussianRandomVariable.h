#ifndef _GaussianRandomVariable_h
#define _GaussianRandomVariable_h


#include "AbstractRandomVariable.h"

class GaussianRandomVariable : public AbstractRandomVariable{
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    GaussianRandomVariable(double mean, double variance);
    ~GaussianRandomVariable();

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Getter(s) and Setter(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    virtual ScalarType GetParameter(std::string param_name) const;
    virtual ScalarType GetParameter(int param_key) const;

    inline double GetVariance() const { return variance_; }

    inline double GetMean() const { return mean_; }

    inline void SetMean(double mean) { mean_ = mean; };

    inline void SetVariance(double variance) {variance_ = variance; };

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Draw a sample
    virtual double Sample();

    /// Compute the likelihood given a current state
    virtual double Likelihood(double x);

    /// Compute the loglikelihood given a current state
    virtual double LogLikelihood(double x);

    /// Update the random variable parameters
    virtual void Update(StringScalarHash params);
    virtual void Update(IntScalarHash    params);

protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// mean of the Gaussian
    double mean_;

    /// variance of the Gaussian
    double variance_;

};


#endif //_GaussianRandomVariable_h
