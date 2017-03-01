#include "UniformRandomVariable.h"

#include "GaussianRandomVariable.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

UniformRandomVariable
::UniformRandomVariable(double Min, double Max)
{
  // TODO : Verify '>' or '>=' but I guess a "Dirac Uniform distribution" is unknown
  assert(m_Max > m_Min);
  m_Min = Min;
  m_Max = Max;
}

UniformRandomVariable
::~UniformRandomVariable()
{ }


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Getter(s)  and Setter(s):
////////////////////////////////////////////////////////////////////////////////////////////////////

ScalarType 
UniformRandomVariable
::GetParameter(std::string ParameterName) const 
{
  if(ParameterName == "Min")
      return m_Min;
  else if(ParameterName == "Max")
      return m_Max;
  else
      std::cerr << "This Parameter does not exist";
}


ScalarType 
UniformRandomVariable
::GetParameter(int ParameterKey) const 
{
  if(ParameterKey == 0)
      return m_Min;
  else if(ParameterKey == 1)
      return m_Max;
  else
      std::cerr << "This Parameter does not exist";
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

double
UniformRandomVariable
::Sample()
{
  // TODO : Add in the attribute
  std::uniform_real_distribution<double> Distribution(m_Min, m_Max);

  double Sample =  Distribution(Generator);
  return Sample;
}

double
UniformRandomVariable
::Likelihood(double X)
{
  return 1./(m_Max - m_Min);
}


double 
UniformRandomVariable
::LogLikelihood(double X) 
{
  return log(1./(m_Max - m_Min));
  
}

void
UniformRandomVariable
::Update(StringScalarHash Parameters) 
{
  bool FindAnything = false;
  
  if(Parameters.find("Min") != Parameters.end())
  {
      m_Min = Parameters.at("Min");
      FindAnything = true;
  }
  if(Parameters.find("Max") != Parameters.end())
  {
      m_Max = Parameters.at("Max");
      FindAnything = true;
  }
  
  assert(m_Max > m_Min);
  
  if(!FindAnything) {std::cerr << "The random variable parameter to update does not exist"; }
}

void
UniformRandomVariable
::Update(IntScalarHash Parameters) 
{
  bool FindAnything = false;
  
  if(Parameters.find(0) != Parameters.end())
  {
      m_Min = Parameters.at(0);
      FindAnything = true;
  }
  if(Parameters.find(1) != Parameters.end())
  {
      m_Max = Parameters.at(1);
      FindAnything = true;
  }
  
  assert(m_Max > m_Min);
  
  if(!FindAnything) {std::cerr << "The random variable parameter to update does not exist"; }
}