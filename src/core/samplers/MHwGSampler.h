#pragma once

#include "AbstractSampler.h"

class MHwGSampler : public AbstractSampler {
public:
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  MHwGSampler();
  virtual ~MHwGSampler();
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initialize the sampler
  virtual void InitializeSampler(std::shared_ptr<CandidateRandomVariables>& candidates, AbstractModel& model);

  /// Sample new realizations
  virtual void Sample(Realizations& reals, AbstractModel& model, const Observations& obs);

  
 private:
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  MHwGSampler(const MHwGSampler&);
  MHwGSampler& operator=(const MHwGSampler&);
  
  /// Sample one realization
  void MonoSample(int block_num, Realizations& reals, AbstractModel& model, const Observations& obs);
  
  /// Compute Prior part of the ratio while updating the realization
  ScalarType ComputePriorRatioAndUpdateRealizations(Realizations& reals, const AbstractModel& model);
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /// Uniform distrubution to compute the acceptance rate
  std::uniform_real_distribution<double> uniform_distrib_;
  
  /// Current iteration Number
  int cur_iter_;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Attribute related to the on-going block sampling
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /// Current Block Info
  MiniBlock cur_block_info_;
  
  /// Current realization name
  std::string cur_real_name_;
  
  /// Current realization key
  int cur_real_key_;
  
  /// Current realization number
  int cur_real_num_;
  
  /// Current realization value
  ScalarType cur_real_value_;
  

  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Adaptive-sampling related attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////
    
  /// Acceptance ratios
  IntVectorHash acceptance_rates_;
  
  /// Expected acceptance rate goal
  ScalarType acceptance_rate_goal;
  
  /// Bandwitdh / window to compute the acceptance rate on 
  unsigned int acceptance_rate_window_;
  
  /// Time step when the acceptance rate is computed, in order to increase/decrease the variance of the proposition distribution  
  unsigned int acceptance_rate_computation_;
  
  
};
