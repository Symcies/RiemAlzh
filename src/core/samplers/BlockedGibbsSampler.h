#pragma once

#include <algorithm>

#include "AbstractSampler.h"

class BlockedGibbsSampler : public AbstractSampler {

public:
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  BlockedGibbsSampler();
  virtual ~BlockedGibbsSampler();
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initialize the sampler
  virtual void InitializeSampler(std::shared_ptr<CandidateRandomVariables>& candidates, AbstractModel& model);

  /// Sample new realizations
  virtual void Sample(Realizations& reals, AbstractModel& models, const Observations& obs);


private:
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  BlockedGibbsSampler(const BlockedGibbsSampler&);
  BlockedGibbsSampler& operator=(const BlockedGibbsSampler&);

  /// Sample one block
  void OneBlockSample(int block_num, Realizations& reals, AbstractModel& models, const Observations& obs);

  /// Compute Prior part of the ratio while updating the realization
  ScalarType ComputePriorRatioAndUpdateRealizations(Realizations& reals, const AbstractModel& models);
  
  /// Compute likelihood based on the block type
  VectorType ComputeLogLikelihood(AbstractModel& models, const Observations& obs);

  /// Get previously computed log likelihood
  double GetPreviousLogLikelihood(AbstractModel& model);

  /// Update the last log likelihood computed
  void UpdateLastLogLikelihood(AbstractModel& model, VectorType& computed_log_likelihood);



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Current Block Info
  MiniBlock cur_block_info_;

  /// Current iteration Number
  int cur_iter_;

  /// Parameters updated in the block
  std::vector<std::string> cur_block_params_;

  /// Parameters to recover back in the realization
  std::vector<std::tuple<std::string, unsigned int, ScalarType >> recover_params_;

  /// Uniform distrubution
  std::uniform_real_distribution<double> uniform_distrib_;

};

