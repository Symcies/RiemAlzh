#ifndef _GaussianRandomvariable_h
#define _GaussianRandomvariable_h


#include "AbstractRandomVariable.h"


class GaussianRandomVariable : public AbstractRandomVariable{
public:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    GaussianRandomVariable(double Mean, double Variance);


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Getter(s) and Setter(s):
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Get the mean of the gaussian distribution
    inline double GetMean() const {return m_Mean;};

    // Get the variance of the gaussian distribution
    inline double GetVariance() const { return m_Variance; };

    // Set the mean of the gaussian distribution
    inline void SetMean(double Mean) {m_Mean = Mean;};

    // Set the variance of the gaussian distribution
    inline void SetVariance(double Variance) { m_Variance = Variance; };

    // Get the density estimation
    inline double GetDensity() const;



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Sample to update the current state
    void Sample();


protected:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Mean of the Gaussian
    double m_Mean;

    // Variance of the Gaussian
    double m_Variance;

};


#endif //_GaussianRandomvariable_h
