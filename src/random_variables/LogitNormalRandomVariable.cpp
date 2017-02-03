#include "LogitNormalRandomVariable.h"


#include "GaussianRandomVariable.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

LogitNormalRandomVariable
::LogitNormalRandomVariable(double Mean, double Variance)
{
    m_Mean = Mean;
    m_Variance = Variance;
}

LogitNormalRandomVariable
::~LogitNormalRandomVariable()
{ }


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

double
LogitNormalRandomVariable
::Sample()
{
    std::normal_distribution<double> Distribution(m_Mean, sqrt(m_Variance));
    double Sample = Distribution(Generator);

    return 1.0/ ( 1.0 + std::exp(-Sample));

}

double
LogitNormalRandomVariable
::Likelihood(double X)
{
    if(X < 0 || X > 1) return 0;
    double denom =   sqrt(2.0*m_Variance*M_PI) * X * (1-X);
    double LogitX = std::log( X / (1.0 - X));
    double num = exp( - (LogitX - m_Mean)*(LogitX - m_Mean) / (2.0*m_Variance));
    double result = num/denom;


    if(isnan(result)) std::cout << "Mean/Variance : " << m_Mean << "/" << m_Variance << "  -  Num/Denom : "<<  num << "/" << denom << std::endl;
    return result;
}


double 
LogitNormalRandomVariable
::LogLikelihood(double X) 
{
    double LogLikelihood = std::log(Likelihood(X));
    return LogLikelihood;
    
}