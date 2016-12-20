#include "LogNormalRandomVariable.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////
    

LogNormalRandomVariable
::LogNormalRandomVariable(double Mean, double Variance) 
{
    m_Mean = Mean;
    m_Variance = Variance;
}

LogNormalRandomVariable
::~LogNormalRandomVariable() 
{}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

double 
LogNormalRandomVariable
::Sample() 
{
    std::random_device RD;
    std::mt19937_64 Generator(RD());
    std::lognormal_distribution<double> Distribution(m_Mean, sqrt(m_Variance));
    
    return Distribution(Generator);
}


double 
LogNormalRandomVariable
::Likelihood(double X) 
{
    double denom = sqrt(2.0*m_Variance*M_PI) * X;
    double num = exp( - (std::log(X) - m_Mean) * (std::log(X) - m_Mean) / (2.0 * m_Variance ));
    double result = num/denom;
    
    if(isnan(result)) std::cout << "Mean/Variance : " << m_Mean << "/" << m_Variance << " - Num/Denom : " << num << "/" << denom << std::endl;
    return result;
}


double 
LogNormalRandomVariable
::LogLikelihood(double X) 
{
    double LogLikelihood = - 1.0/2.0 * log(2.0*m_Variance*M_PI) - std::log(X);
    LogLikelihood +=  - (std::log(X) - m_Mean)*(std::log(X) - m_Mean) / (2.0 * m_Variance);
    return LogLikelihood;   
}



