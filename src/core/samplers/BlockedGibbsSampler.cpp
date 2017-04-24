#include "BlockedGibbsSampler.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

BlockedGibbsSampler::BlockedGibbsSampler()
{
  cur_iter_ = 0;
  cur_block_info_ = {std::make_tuple<int, std::string, int>(-1, "All", 0)};
  uniform_distrib_ = std::uniform_real_distribution<double>(0.0, 1.0);
}


BlockedGibbsSampler::~BlockedGibbsSampler()
{

}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void BlockedGibbsSampler::InitializeSampler(std::shared_ptr<CandidateRandomVariables>& candidates, AbstractModel& model)
{
  blocks_ = model.GetSamplerBlocks();
  candidate_rand_var_ = candidates;
}


void BlockedGibbsSampler::Sample(Realizations& reals, AbstractModel& model, const Observations& obs)
{
  /// Initialization of the sampler and the model
  cur_block_info_ = {std::make_tuple<int, std::string, int>(-1, "All", 0)};
  model.UpdateModel(reals, cur_block_info_);
  //std::cout << " 1. " << std::endl;
  VectorType log_likelihood = ComputeLogLikelihood(model, obs);
  //std::cout << "ok --> " << log_likelihood.sum() << std::endl;
  UpdateLastLogLikelihood(model, log_likelihood);

  for (int i = 0; i < blocks_.size(); ++i) {
    OneBlockSample(i, reals, model, obs);
  }

  ++cur_iter_;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void BlockedGibbsSampler::OneBlockSample(int block_num, Realizations& reals,
  AbstractModel &model, const Observations& obs)
{
  /// Initialization
  cur_block_info_ = blocks_.at(block_num);
  cur_block_params_.clear();
  recover_params_.clear();


  /// Loop over the realizations of the block to update the ratio and the realizations
  double acceptation_ratio = ComputePriorRatioAndUpdateRealizations(reals, model);

  /// Compute the previous log likelihood
  acceptation_ratio -= GetPreviousLogLikelihood(model);

  /// Compute the candidate log likelihood
  model.UpdateModel(reals, cur_block_info_, cur_block_params_);
  VectorType computed_log_likelihood = ComputeLogLikelihood(model, obs);
  acceptation_ratio += computed_log_likelihood.sum();

  
  // TODO : TO ERASE
  //std::cout << "   -  LL : " << computed_log_likelihood.sum() << " - " << GetPreviousLogLikelihood(model)  << std::endl;
  
  /// Compute the aceceptance ratio
  acceptation_ratio = exp(std::min(acceptation_ratio, 0.0));

  ///  Rejection : Candidate not accepted
  if (uniform_distrib_(generator) > acceptation_ratio) {
    for(auto it = recover_params_.begin(); it != recover_params_.end(); ++it)
    {
      std::string name = std::get<0>(*it);
      unsigned int real_num = std::get<1>(*it);
      ScalarType prev_real = std::get<2>(*it);
      reals.at(name, real_num) = prev_real;
    }

    //std::cout << "     -  REJECT -  likelihood : " << GetPreviousLogLikelihood(model) << std::endl; 
    model.UpdateModel(reals, cur_block_info_, cur_block_params_);
  }
      /// Acceptation : Candidate is accepted
  else {
    //std::cout << "     -  ACCEPT -  likelihood : " << computed_log_likelihood.sum() << std::endl; 
    UpdateLastLogLikelihood(model, computed_log_likelihood);
  }
  
}

ScalarType BlockedGibbsSampler::ComputePriorRatioAndUpdateRealizations(Realizations& reals, const AbstractModel& model)
{
  double acceptance_ratio = 0;

  for (auto it = cur_block_info_.begin(); it != cur_block_info_.end(); ++it) {
    /// Initialization
    std::string name_real = std::get<1>(*it);
    int key = reals.ReverseNameToKey(name_real);
    unsigned int real_num = std::get<2>(*it);
    cur_block_params_.push_back(name_real);

    /// Get the current realization and recover it
    auto cur_rand_var = model.GetRandomVariable(key);
    ScalarType cur_real = reals.at(key, real_num);
    recover_params_.push_back(std::make_tuple(name_real, real_num, cur_real));

    /// Get a candidate realization
    auto candidate_rand_var = candidate_rand_var_->GetRandomVariable(key, real_num);
    candidate_rand_var.SetMean(cur_real);
    ScalarType candidate_real = candidate_rand_var.Sample();

    /// Calculate the acceptance ratio
    /*
    acceptance_ratio += cur_rand_var->LogLikelihood(candidate_real);
    acceptance_ratio -= cur_rand_var->LogLikelihood(cur_real);
    */
    ScalarType ok1 = cur_rand_var->LogLikelihood(candidate_real);
    ScalarType ok2 = cur_rand_var->LogLikelihood(cur_real);
    acceptance_ratio += ok1 - ok2;
    //std::cout << "Prior : " << ok1 << " - " << ok2;
     
    /// Update the NewRealizations
    reals.at(name_real, real_num) = candidate_real;
    
    
    //std::cout << "if " << cur_real - candidate_real << " > 0, then  "; 
  }

  return acceptance_ratio;

}



BlockedGibbsSampler::VectorType BlockedGibbsSampler::ComputeLogLikelihood(AbstractModel& model, const Observations& obs)
{
  return model.ComputeLogLikelihood(obs, cur_block_info_);
}

double BlockedGibbsSampler::GetPreviousLogLikelihood(AbstractModel& model)
{
  return model.GetPreviousLogLikelihood(cur_block_info_);
}

void BlockedGibbsSampler::UpdateLastLogLikelihood(AbstractModel& model, VectorType& computed_log_likelihood)
{
  model.SetPreviousLogLikelihood(computed_log_likelihood, cur_block_info_);

}
