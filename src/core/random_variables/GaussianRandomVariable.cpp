#include "GaussianRandomVariable.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

GaussianRandomVariable::GaussianRandomVariable(double mean, double variance)
{
  mean_ = mean;
  variance_ = variance;
}

GaussianRandomVariable::GaussianRandomVariable(const GaussianRandomVariable& gr_var)
{
  mean_ = gr_var.mean_;
  variance_ = gr_var.variance_;
}

GaussianRandomVariable::~GaussianRandomVariable()
{ }

GaussianRandomVariable& GaussianRandomVariable::operator=(const GaussianRandomVariable& gr_var)
{
  mean_ = gr_var.mean_;
  variance_ = gr_var.variance_;

  return *this;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Getter(s)  and Setter(s):
////////////////////////////////////////////////////////////////////////////////////////////////////

ScalarType GaussianRandomVariable::GetParameter(std::string param_name) const
{
  if(param_name == "Mean")
      return mean_;
  else if(param_name == "Variance")
      return variance_;
  else
      std::cerr << "This Parameter does not exist";
}


ScalarType GaussianRandomVariable::GetParameter(int param_key) const
{
  if(param_key == 0)
      return mean_;
  else if(param_key == 1)
      return variance_;
  else
      std::cerr << "This Parameter does not exist";
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

double GaussianRandomVariable::Sample()
{
  // TODO : Add in the attribute !
  std::normal_distribution<double> distibution(mean_, sqrt(variance_));

  double sample =  distibution(generator);
  return sample;
}

double GaussianRandomVariable::Likelihood(double x)
{
  double denom =   sqrt(2.0 * variance_ * M_PI);
  double num = exp( - (x - mean_)*(x - mean_) / (2.0*variance_));
  double result = num/denom;


  if(isnan(result)) std::cout << "mean/variance : " << mean_ << "/" << variance_ << "  -  Num/Denom : "<<  num << "/" << denom << std::endl;
  return result;
}


double GaussianRandomVariable::LogLikelihood(double x)
{
  double log_likelihood = - 1.0/2.0 * log(2.0*variance_*M_PI);
  log_likelihood +=  - (x - mean_)*(x - mean_) / (2.0 * variance_);
  return log_likelihood;

}

void GaussianRandomVariable::Update(StringScalarHash params)
{
  bool find_anything = false;

  if(params.find("Mean") != params.end())
  {
      mean_ = params.at("Mean");
      find_anything = true;
  }
  if(params.find("Variance") != params.end())
  {
      assert(params.at("Variance") > 0);
      variance_ = params.at("Variance");
      find_anything = true;
  }

  if(!find_anything) {std::cerr << "The random variable parameter to update does not exist"; }
}

void GaussianRandomVariable::Update(IntScalarHash params)
{
  bool find_anything = false;

  if(params.find(0) != params.end())
  {
      mean_ = params.at(0);
      find_anything = true;
  }
  if(params.find(1) != params.end())
  {
      assert(params.at(1) > 0);
      variance_ = params.at(1);
      find_anything = true;
  }

  if(!find_anything) {std::cerr << "The random variable parameter to update does not exist"; }
}
