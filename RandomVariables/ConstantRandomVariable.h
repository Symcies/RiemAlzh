#ifndef Longitudinal_ConstantRandomVariable_h
#define Longitudinal_ConstantRandomVariable_h

#include <iostream>
#include <cmath>
#include "AbstractRandomVariable.h"

class ConstantRandomVariable : public AbstractRandomVariable {
public:
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    /// Constructor 
    ConstantRandomVariable(double Mean);
    
    /// Copy constructor
    ConstantRandomVariable(const ConstantRandomVariable& CRV);
    
    
    
    /// Destructor
    ~ConstantRandomVariable();

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Getter(s) and Setter(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    inline double GetMean() { return m_Mean; }
    
    inline void SetMean(double Mean) { m_Mean = Mean; }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Draw a new sample
    virtual double Sample();

    /// Compute the likelihood
    virtual double Likelihood(double X);

    /// Compute the log likelihood
    virtual double LogLikelihood(double X);


protected:

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Constant value of the random variable
    double m_Mean;
};


#endif //Longitudinal_ConstantRandomVariable_h
