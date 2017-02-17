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
/// Getter(s)  and Setter(s):
////////////////////////////////////////////////////////////////////////////////////////////////////

ScalarType 
GaussianRandomVariable
::GetParameter(std::string ParameterName) const 
{
    if(ParameterName == "Mean")
        return m_Mean;
    else if(ParameterName == "Variance")
        return m_Variance;
    else
        std::cerr << "This Parameter does not exist";
}


ScalarType 
GaussianRandomVariable
::GetParameter(int ParameterKey) const 
{
    if(ParameterKey == 0)
        return m_Mean;
    else if(ParameterKey == 1)
        return m_Variance;
    else
        std::cerr << "This Parameter does not exist";
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
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

void
GaussianRandomVariable
::Update(StringScalarHash Parameters) 
{
    bool FindAnything = false;
    
    if(Parameters.find("Mean") != Parameters.end())
    {
        m_Mean = Parameters.at("Mean");
        FindAnything = true;
    }
    if(Parameters.find("Variance") != Parameters.end())
    {
        assert(Parameters.at("Variance") > 0);
        m_Variance = Parameters.at("Variance");
        FindAnything = true;
    }
    
    if(!FindAnything) {std::cerr << "The random variable parameter to update does not exist"; }
}

void
GaussianRandomVariable
::Update(IntScalarHash Parameters) 
{
    bool FindAnything = false;
    
    if(Parameters.find(0) != Parameters.end())
    {
        m_Mean = Parameters.at(0);
        FindAnything = true;
    }
    if(Parameters.find(1) != Parameters.end())
    {
        assert(Parameters.at(1) > 0);
        m_Variance = Parameters.at(1);
        FindAnything = true;
    }
    
    if(!FindAnything) {std::cerr << "The random variable parameter to update does not exist"; }
}