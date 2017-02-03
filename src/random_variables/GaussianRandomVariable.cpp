#include "GaussianRandomVariable.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

GaussianRandomVariable
::GaussianRandomVariable(double Mean, double Variance)
{
    m_Mean = Mean;
    m_Variance = Variance;
}

GaussianRandomVariable
::~GaussianRandomVariable()
{ }

////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

double
GaussianRandomVariable
::Sample()
{
    std::normal_distribution<double> Distribution(m_Mean, sqrt(m_Variance));

    double Sample =  Distribution(Generator);
    return Sample;
}

double
GaussianRandomVariable
::Likelihood(double X)
{
    double denom =   sqrt(2.0*m_Variance*M_PI);
    double num = exp( - (X - m_Mean)*(X - m_Mean) / (2.0*m_Variance));
    double result = num/denom;


    if(isnan(result)) std::cout << "Mean/Variance : " << m_Mean << "/" << m_Variance << "  -  Num/Denom : "<<  num << "/" << denom << std::endl;
    return result;
}


double 
GaussianRandomVariable
::LogLikelihood(double X) 
{
    double LogLikelihood = - 1.0/2.0 * log(2.0*m_Variance*M_PI);
    LogLikelihood +=  - (X - m_Mean)*(X - m_Mean) / (2.0 * m_Variance);
    return LogLikelihood;
    
}