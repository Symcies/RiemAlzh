#include "UniformRandomVariable.h"

#include "GaussianRandomVariable.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

UniformRandomVariable::UniformRandomVariable(double min, double max)
{
  // TODO : Verify '>' or '>=' but I guess a "Dirac Uniform distribution" is unknown
  min_ = min;
  max_ = max;
  assert(max_ >= min_);
}

UniformRandomVariable::~UniformRandomVariable()
{ }

UniformRandomVariable::UniformRandomVariable(const UniformRandomVariable& ur_var)
{
  min_ = ur_var.min_;
  max_ = ur_var.max_;
}

UniformRandomVariable& UniformRandomVariable::operator=(const UniformRandomVariable& ur_var)
{
  min_ = ur_var.min_;
  max_ = ur_var.max_;
  return *this;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Getter(s)  and Setter(s):
////////////////////////////////////////////////////////////////////////////////////////////////////

ScalarType UniformRandomVariable::GetParameter(std::string param_name) const
{
  if(param_name == "min")
      return min_;
  else if(param_name == "max")
      return max_;
  else
      std::cerr << "This Parameter does not exist";
}


ScalarType UniformRandomVariable::GetParameter(int param_key) const
{
  if(param_key == 0)
      return min_;
  else if(param_key == 1)
      return max_;
  else
      std::cerr << "This Parameter does not exist";
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

double UniformRandomVariable::Sample()
{
  // TODO : Add in the attribute
  std::uniform_real_distribution<double> distibution(min_, max_);

  double sample =  distibution(generator);
  return sample;
}

double UniformRandomVariable::Likelihood(double x)
{
  return 1./(max_ - min_);
}


double UniformRandomVariable::LogLikelihood(double x)
{
  return log(1./(max_ - min_));

}

void UniformRandomVariable::Update(StringScalarHash params)
{
  bool find_anything = false;

  if(params.find("min") != params.end())
  {
      min_ = params.at("min");
      find_anything = true;
  }
  if(params.find("max") != params.end())
  {
      max_ = params.at("max");
      find_anything = true;
  }

  assert(max_ > min_);

  if(!find_anything) {std::cerr << "The random variable parameter to update does not exist"; }
}

void UniformRandomVariable::Update(IntScalarHash params)
{
  bool find_anything = false;

  if(params.find(0) != params.end())
  {
      min_ = params.at(0);
      find_anything = true;
  }
  if(params.find(1) != params.end())
  {
      max_ = params.at(1);
      find_anything = true;
  }

  assert(max_ > min_);

  if(!find_anything) {std::cerr << "The random variable parameter to update does not exist"; }
}
