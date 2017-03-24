#include "BlockedGibbsSampler.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


BlockedGibbsSampler::BlockedGibbsSampler()
{
  cur_iter_ = 0;
  cur_block_type_ = -1;
  uniform_distrib_ = std::uniform_real_distribution<double>(0.0, 1.0);
}

BlockedGibbsSampler::BlockedGibbsSampler(unsigned int memoryless_sampling_time, double expected_acceptance_ratio)
{
  uniform_distrib_ = std::uniform_real_distribution<double>(0.0, 1.0);
  memoryless_sampling_time_ = memoryless_sampling_time;
  expected_acceptance_ratio_ = expected_acceptance_ratio;
  cur_iter_ = 0;
  cur_block_type_ = -1;
}

BlockedGibbsSampler::~BlockedGibbsSampler()
{

}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void BlockedGibbsSampler::InitializeSampler(Realizations& reals, AbstractModel &model)
{
  candidate_rand_var_.InitializeCandidateRandomVariables(reals, model);
  cur_block_type_ = -1;
  blocks_ = model.GetSamplerBlocks();

}


void BlockedGibbsSampler::Sample(Realizations& reals, AbstractModel& model, const Observations& obs)
{
  cur_block_type_ = -1;
  ////////////////////////////////////////
  // TODO : Check if the update is needed -> Yes, needed because the orthonormal basis,
  // or the A matrix may have changed because some variables such as V0 or P0 have changed
  model.UpdateModel(reals, -1);
  VectorType log_likelihood = ComputeLogLikelihood(model, obs);
  UpdateLastLogLikelihood(model, log_likelihood);
  ////////////////////////////////////////


  for (int i = 0; i < blocks_.size(); ++i)
  {
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
  SamplerBlock current_block = blocks_.at(block_num);
  cur_block_type_ = current_block.first;
  cur_block_params_.clear();
  cur_block_params_bis_.clear();
  recover_params_.clear();

  /// Loop over the realizations of the block to update the ratio and the realizations
  double acceptation_ratio = ComputePriorRatioAndUpdateRealizations(reals, model, current_block.second);

  /// Compute the previous log likelihood
  acceptation_ratio -= GetPreviousLogLikelihood(model);

  /// Compute the candidate log likelihood
  model.UpdateModel(reals, cur_block_type_, cur_block_params_);
  VectorType computed_log_likelihood = ComputeLogLikelihood(model, obs);
  acceptation_ratio += computed_log_likelihood.sum();

  /// Compute the aceceptance ratio
  acceptation_ratio = std::min(acceptation_ratio, 0.0);
  acceptation_ratio = exp(acceptation_ratio);

  ///  Rejection : Candidate not accepted
  if(uniform_distrib_(generator) > acceptation_ratio)
  {
    for(auto it = recover_params_.begin(); it != recover_params_.end(); ++it)
    {
      std::string name = std::get<0>(*it);
      unsigned int real_num = std::get<1>(*it);
      ScalarType prev_real = std::get<2>(*it);
      reals.at(name, real_num) = prev_real;
    }

    model.UpdateModel(reals, cur_block_type_, cur_block_params_);
  }
      /// Acceptation : Candidate is accepted
  else
  {
    UpdateLastLogLikelihood(model, computed_log_likelihood);
  }

  /// Adaptative variances for the realizations
  UpdateBlockRandomVariable(acceptation_ratio, current_block .second);

}

ScalarType BlockedGibbsSampler::ComputePriorRatioAndUpdateRealizations(Realizations& reals, const AbstractModel& model, const MiniBlock& var)
{
  double acceptance_ratio = 0;

  for(auto it = var.begin(); it != var.end(); ++it)
  {
      /// Initialization
      std::string name_real = it->first;
      int key = reals.ReverseNameToKey(name_real);
      unsigned int real_num = it->second;
      cur_block_params_.push_back(name_real);
      cur_block_params_bis_.push_back(key);

      /// Get the current realization and recover it
      auto cur_rand_var = model.GetRandomVariable(key);
      ScalarType cur_real = reals.at(key, real_num);
      //recover_params_[name_real] = {real_num, cur_real};
      recover_params_.push_back(std::make_tuple(name_real, real_num, cur_real));

      /// Get a candidate realization
      auto candidate_rand_var = candidate_rand_var_.GetRandomVariable(key, real_num);
      candidate_rand_var.SetMean(cur_real);
      ScalarType candidate_real = candidate_rand_var.Sample();

      /// Calculate the acceptance ratio
      acceptance_ratio += cur_rand_var->LogLikelihood(candidate_real);
      acceptance_ratio -= cur_rand_var->LogLikelihood(cur_real);

      /// Update the NewRealizations
      reals.at(name_real, real_num) = candidate_real;

  }

  return acceptance_ratio;

}



BlockedGibbsSampler::VectorType BlockedGibbsSampler::ComputeLogLikelihood(AbstractModel& model, const Observations& obs)
{
  return model.ComputeLogLikelihood(obs, cur_block_type_);
}

double BlockedGibbsSampler::GetPreviousLogLikelihood(AbstractModel& model)
{
  return model.GetPreviousLogLikelihood(cur_block_type_);
}

void BlockedGibbsSampler::UpdateLastLogLikelihood(AbstractModel& model, VectorType& computed_log_likelihood)
{
  model.SetPreviousLogLikelihood(computed_log_likelihood, cur_block_type_);

}


void BlockedGibbsSampler::UpdateBlockRandomVariable(double acceptance_ratio, const MiniBlock &var)
{
  double epsilon = DecreasingStepSize(cur_iter_, memoryless_sampling_time_);
  double denom = expected_acceptance_ratio_;

  if(acceptance_ratio > expected_acceptance_ratio_)
      denom = 1 - expected_acceptance_ratio_;

  for(auto it = var.begin(); it != var.end(); ++it)
  {
      /// Initialize
      std::string name_real = it->first;
      unsigned int real_num = it->second;

      /// Compute newvariance
      GaussianRandomVariable gaussian_rand_var = candidate_rand_var_.GetRandomVariable(name_real, real_num);
      double cur_variance = gaussian_rand_var.GetVariance();
      double new_variance = cur_variance * (1 + epsilon * (acceptance_ratio - expected_acceptance_ratio_) / denom );
      new_variance = std::max(new_variance, 0.0000000000000000000000000000001);
      new_variance = std::max(new_variance, cur_variance / 50);
      new_variance = std::min(new_variance, cur_variance * 50);

      /// Set new variance
      //candidate_rand_var_.UpdatePropositionVariableVariance(name_real, real_num, new_variance);

  }
}
