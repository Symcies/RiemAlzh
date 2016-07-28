#ifndef _GaussianRandomVariable_h
#define _GaussianRandomVariable_h


#include "AbstractRandomVariable.h"

class GaussianRandomVariable : public AbstractRandomVariable{
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    GaussianRandomVariable(double Mean, double Variance);
    ~GaussianRandomVariable();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    inline double GetVariance() { return m_Variance; }

    inline void SetMean(double Mean) { m_Mean = Mean; };

    inline void SetVariance(double Variance) {m_Variance = Variance; };


    /// Draw a sample
    virtual double Sample();

    /// Compute the likelihood given a current state
    virtual double Likelihood(double X);



protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Mean of the Gaussian
    double m_Mean;

    /// Variance of the Gaussian
    double m_Variance;

};


#endif //_GaussianRandomVariable_h
