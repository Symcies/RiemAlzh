#include "GaussianRandomVariable.h"


GaussianRandomVariable
::GaussianRandomVariable(double Mean, double Variance)
{
    m_Mean = Mean;
    m_Variance = Variance;
    Sample();
}

void
GaussianRandomVariable
::Sample()
{
    std::random_device RD;
    std::mt19937_64 Generator(RD());
    std::normal_distribution<double> Distribution(m_Mean, m_Variance);

    m_CurrentState = Distribution(Generator);
    //std::cout << " Mean / Variance / Current : " << m_Mean << "  /  " << m_Variance << "  /  " << m_CurrentState << std::endl;
}

double
GaussianRandomVariable
::GetDensity()
const
{
    double denom =   sqrt(2*m_Variance*M_PI);
    double num = exp( - (m_CurrentState - m_Mean)*(m_CurrentState - m_Mean) / (2*m_Variance));
    double result = num/denom;
    if(isnan(result)) std::cout << "Mean/Variance : " << m_Mean << "/" << m_Variance << "  -  Num/Denom : "<<  num << "/" << denom << std::endl;
    return result;
}

