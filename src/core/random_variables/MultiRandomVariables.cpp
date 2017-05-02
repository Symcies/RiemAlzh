#include "MultiRandomVariables.h"




////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

MultiRandomVariables::MultiRandomVariables()
{

}

MultiRandomVariables::~MultiRandomVariables()
{

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<AbstractRandomVariable> MultiRandomVariables::GetRandomVariable(std::string name) const
{
  if(rand_var_string_type_key_.at(string_to_int_key_.at(name)) == "Gaussian")
  {
    // TODO : Change to GetParameters(0) and GetPArameters(1)
    ScalarType mean = rand_var_.at(string_to_int_key_.at(name))->GetParameter("Mean");
    ScalarType variance = rand_var_.at(string_to_int_key_.at(name))->GetParameter("Variance");
    return std::make_unique<GaussianRandomVariable>(mean, variance);
  }
}

std::unique_ptr<AbstractRandomVariable> MultiRandomVariables::GetRandomVariable(int key)
const
{
  if(rand_var_int_type_key_.at(key) == 0)
  {
    // TODO : Change to GetParameters(0) and GetPArameters(1)
    ScalarType mean = rand_var_.at(key)->GetParameter("Mean");
    ScalarType variance = rand_var_.at(key)->GetParameter("Variance");
    return std::make_unique<GaussianRandomVariable>(mean, variance);
  }
}

void MultiRandomVariables::Clear()
{
  rand_var_.clear();
  string_to_int_key_.clear();
  int_to_string_key_.clear();
  key_count_ = 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void MultiRandomVariables::AddRandomVariable(std::string name, std::string type, const std::vector<double>& params)
{
  ////////////////////////////////////////////////////////////
  //// The rand_var is not in the random variables
  ////////////////////////////////////////////////////////////
  if(string_to_int_key_.find(name) == string_to_int_key_.end())
  {
    //std::cout << "You add the random variable " << name << std::endl;
    std::shared_ptr<AbstractRandomVariable> rand_var;
    if (type == "Gaussian") {
      rand_var = std::make_shared<GaussianRandomVariable>(params[0], params[1]);
      rand_var_string_type_key_.insert({key_count_, "Gaussian"});
      rand_var_int_type_key_.insert({key_count_, 0});
    }

    string_to_int_key_.insert({name, key_count_});
    int_to_string_key_.insert({key_count_, name});
    rand_var_.insert({key_count_, rand_var});
    ++key_count_;
  }
  ////////////////////////////////////////////////////////////
  //// The rand_var is already in the random variables of the class
  ////////////////////////////////////////////////////////////
  else
  {
    int key = string_to_int_key_.at(name);
    std::string prev_type = rand_var_string_type_key_.at(key);
    std::shared_ptr<AbstractRandomVariable> rand_var;

    if(type == prev_type)
    {
      //std::cout << "You overwrite the random variable " << name << " with the same type" << std::endl;
      rand_var = std::make_shared<GaussianRandomVariable>(params[0], params[1]);
    }
    else
    {
      std::cout << "You overwrite the random variable " << name << " with a new type" << std::endl;
      std::cerr << "TODO in the multi random variables ..." << std::endl;
    }

    rand_var_[key] = rand_var;
  }
}

void MultiRandomVariables::UpdateRandomVariable(std::string name, IntScalarHash params)
{
  if(rand_var_.find(string_to_int_key_.at(name)) != rand_var_.end())
    rand_var_.at(string_to_int_key_.at(name))->Update(params);
  else
    std::cerr << name << " is not found in the random variables";
}

void MultiRandomVariables::UpdateRandomVariable(std::string name, StringScalarHash params)
{
  if(rand_var_.find(string_to_int_key_.at(name)) != rand_var_.end())
    rand_var_.at(string_to_int_key_.at(name))->Update(params);
  else
    std::cerr << name << " is not found in the random variables";
}


void MultiRandomVariables::UpdateRandomVariable(int key, IntScalarHash params)
{
  if(rand_var_.find(key) != rand_var_.end())
    rand_var_.at(key)->Update(params);
  else
    std::cerr << int_to_string_key_.at(key) << " is not found in the random variable";
}


void MultiRandomVariables::UpdateRandomVariable(int key, StringScalarHash params)
{
  if(rand_var_.find(key) != rand_var_.end())
    rand_var_.at(key)->Update(params);
  else
    std::cerr << int_to_string_key_.at(key) << " is not found in the random variable";
}


Realizations MultiRandomVariables::SimulateRealizations(StringIntHash num_of_real_per_rand_var)
{
  Realizations realizations;


  for(auto it = num_of_real_per_rand_var.begin(); it != num_of_real_per_rand_var.end(); ++it)
  {
    std::string name = it->first;
    int key = string_to_int_key_.at(name);
    int reals_num = it->second;

    VectorType real_vec = rand_var_.at(key)->Samples(reals_num);
    realizations.AddRealizations(name, key, real_vec);
  }


  return realizations;
}
