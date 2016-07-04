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
::SetMean(double Mean)
{
    m_Mean = Mean;
}

void
GaussianRandomVariable
::SetVariance(double Variance)
{
    m_Variance = Variance;
}

void
GaussianRandomVariable
::Sample()
{
    std::random_device RD;
    std::default_random_engine Generator(RD());
    std::normal_distribution<double> Distribution(m_Mean, m_Variance);

    m_CurrentState = Distribution(Generator);
}

double
GaussianRandomVariable
::GetDensity()
const
{
    double denom =   sqrt(2*m_Variance*M_PI);
    double num = exp( - (m_CurrentState - m_Mean) / (2*m_Variance));
    return num/denom;
}

