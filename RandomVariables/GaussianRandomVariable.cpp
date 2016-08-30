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
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

double
GaussianRandomVariable
::Sample()
{
    std::random_device RD;
    std::mt19937_64 Generator(RD());
    std::normal_distribution<double> Distribution(m_Mean, m_Variance);

    return Distribution(Generator);

}

double
GaussianRandomVariable
::Likelihood(double X)
{
    double denom =   sqrt(2*m_Variance*M_PI);
    double num = exp( - (X - m_Mean)*(X - m_Mean) / (2*m_Variance));
    double result = num/denom;

    //std::cout << "Mean/Variance : " << m_Mean << "/" << m_Variance << std::endl;

    if(isnan(result)) std::cout << "Mean/Variance : " << m_Mean << "/" << m_Variance << "  -  Num/Denom : "<<  num << "/" << denom << std::endl;
    return result;
}