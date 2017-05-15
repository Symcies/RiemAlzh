#include "MHwGSampler.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


MHwGSampler::MHwGSampler() {
  
  cur_iter_ = 0;
  cur_block_info_ = {std::make_tuple<int, std::string, int>(-1, "All", 0)};
  uniform_distrib_ = std::uniform_real_distribution<double>(0.0, 1.0);
  
  // TODO TODO TODO TODO TODO TODO 
  // Put the following arguments in the xml file
  acceptance_rate_computation_ = 10;
  acceptance_rate_window_ = 10;
  assert(acceptance_rate_window_ <= acceptance_rate_computation_);
  acceptance_rate_goal = 0.35;
  
  
  
}


MHwGSampler::~MHwGSampler() {
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void MHwGSampler::InitializeSampler(std::shared_ptr<CandidateRandomVariables> &candidates, AbstractModel& model) {
  candidate_rand_var_ = candidates;
  blocks_ = model.GetSamplerBlocks(1);
  
  // Check if the blocks have only one variable !
  for (auto it = blocks_.begin(); it != blocks_.end(); ++it) {
    assert(it->size() == 1);
  }
   
}

void MHwGSampler::Sample(Realizations &reals, AbstractModel &model, const Observations &obs) {
  /// Initialization of the sampler and the model
  cur_block_info_ = {std::make_tuple<int, std::string, int>(-1, "All", 0)};
  model.UpdateModel(reals, cur_block_info_);
  model.InitializeLogLikelihood(obs);
  
  for (int i = 0; i < blocks_.size(); ++i) {
    MonoSample(i, reals, model, obs);
  }
  
  if(cur_iter_%acceptance_rate_computation_) {
    // TODO : Compute the variance of the proposal distribution
  }

  ++cur_iter_;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Attribute(s)
////////////////////////////////////////////////////////////////////////////////////////////////////

void MHwGSampler::MonoSample(int block_num,
                             Realizations &reals,
                             AbstractModel &model,
                             const Observations &obs) {

  /// Initialize
  cur_block_info_ = blocks_.at(block_num);
  cur_real_name_  = std::get<1>(cur_block_info_[0]);
  cur_real_num_   = std::get<2>(cur_block_info_[0]);
  cur_real_key_   = reals.ReverseNameToKey(cur_real_name_);
  cur_real_value_ = reals.at(cur_real_key_, cur_real_num_);

  /// Compute the prior part of the acceptance rate, and update the realizations
  ScalarType acceptance_rate = ComputePriorRatioAndUpdateRealizations(reals, model);
  
  /// Compute the previous log lokelihood
  acceptance_rate -= model.GetPreviousLogLikelihood(cur_block_info_);
  
  /// Compute the candidate log likelihood
  model.UpdateModel(reals, cur_block_info_, {cur_real_name_});
  VectorType computed_log_likelihood = model.ComputeLogLikelihood(obs, cur_block_info_);
  acceptance_rate += computed_log_likelihood.sum();
  
  /// Compute the aceceptance ratio
  acceptance_rate = exp(std::min(acceptance_rate, 0.0));
  
  ///  Rejection : Candidate not accepted
  if (uniform_distrib_(generator) > acceptance_rate) {
    
    reals.at(cur_real_key_, cur_real_num_) = cur_real_value_;
    model.UpdateModel(reals, cur_block_info_, {cur_real_name_});
    
  } else {
    
    model.SetPreviousLogLikelihood(computed_log_likelihood, cur_block_info_);
    
  }
}


ScalarType MHwGSampler::ComputePriorRatioAndUpdateRealizations(Realizations &reals,
                                                         const AbstractModel &model) {
  
  ScalarType acceptance_rate = 0;
  
  /// Get the current realization and recover it
  auto cur_rand_var   = model.GetRandomVariable(cur_real_key_);
  ScalarType cur_real = reals.at(cur_real_key_, cur_real_num_);
  
  /// Get a candidate realization
  auto candidate_rand_var = candidate_rand_var_->GetRandomVariable(cur_real_key_, cur_real_num_);
  candidate_rand_var.SetMean(cur_real);
  ScalarType candidate_real = candidate_rand_var.Sample();
  
  /// Compute the acceptance rate
  acceptance_rate += cur_rand_var->LogLikelihood(candidate_real);
  acceptance_rate -= cur_rand_var->LogLikelihood(cur_real);
  
  /// Update the realizations
  reals.at(cur_real_key_, cur_real_num_) = candidate_real;
      
  return acceptance_rate;
}