#include "AbstractModel.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


std::shared_ptr< AbstractRandomVariable > AbstractModel::GetRandomVariable(std::string name) const
{
  return rand_var_.GetRandomVariable(name);
}

std::shared_ptr< AbstractRandomVariable > AbstractModel::GetRandomVariable(int key) const
{
  return rand_var_.GetRandomVariable(key);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

ScalarType AbstractModel::InitializePropositionDistributionVariance(std::string name) const {
  name = name.substr(0, name.find_first_of("#"));
  return proposition_distribution_variance_.at(name);
}

Realizations AbstractModel::SimulateRealizations()
{
  return rand_var_.SimulateRealizations(asso_num_real_per_rand_var_);
}
