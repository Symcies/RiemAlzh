#pragma once

#include <algorithm>
#include <memory>
#include <cassert>
#include <cmath>
#include <iostream>

#include "AbstractModel.h"
#include "SamplerSettings.h"
#include "SamplerFactory.h"
#include "AlgorithmSettings.h"

#include "Observations.h"

class Algorithm {
public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  typedef typename LinearAlgebra<ScalarType>::MatrixType MatrixType;
  typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
  typedef typename std::unordered_map<int, VectorType> IntVectorHash;
  typedef std::vector<VectorType> SufficientStatisticsVector;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  Algorithm(io::AlgorithmSettings& settings);
  ~Algorithm();

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  inline void SetModel(std::shared_ptr<AbstractModel>& m ) { model_ = m; }

  inline int GetMaximumNumberOfIterations() { return max_iter_num_; }
  inline int GetNumberOfBurnInIterations()  { return burnin_iter_num_; }
  inline int GetOutputDisplayIteration()    { return output_iter_; }
  inline int GetDataSaveIteration()         { return data_save_iter_; }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Other method :
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /// Compute the MCMC SAEM algorithm
  void ComputeMCMCSAEM(const Observations& obs);
  
  /// Set the samplers 
  void AddSamplers(io::SamplerSettings& settings);
  
private:
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Initialization Method(s):
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /// Initialize the stochastic approximation
  void InitializeStochasticSufficientStatistics(const Observations& obs);

  /// Initialize sampler
  void InitializeSampler();

  /// Initialize Manifold
  void InitializeModel(const Observations& obs);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Method(s) of the MCMC SAEM:
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /// Compute the MCMC SAEM algorithm
  void IterationMCMCSAEM(const Observations& obs, int iter);
  
  /// Compute the simulation step : Gibbs Sampling
  void ComputeSimulationStep(const Observations& obs, int iter);

  /// Compute the stochastic coefficient
  void ComputeStochasticApproximation(SufficientStatisticsVector& stat_vector, int iter);

  /// Compute the decreasing step size of the approximation step
  double DecreasingStepSize(int iter);

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Output Method(s):
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  bool IsOutputIteration(int iter);
  bool IsDataSaveIteration(int iter);
  void DisplayIterations(int iter);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Output(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Number of iterations to wait till next output display
  unsigned int output_iter_;

  /// Number of iterations to wait till next data saving
  unsigned int data_save_iter_;

  /// Acceptance ratios
  IntVectorHash acceptance_ratio_;
  
  /// Random variable whose acceptance ratio are displayed
  std::vector<std::string> acceptance_ratio_to_display_;

  /// Compute the acceptance ratio for each random variable
  void ComputeAcceptanceRatio(Realizations& prev_real_, int iter);

  /// Display acceptance ratio
  void DisplayAcceptanceRatio();

  /// Display Outputs
  void DisplayOutputs();
  
  

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Sampler-related Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /// Get the sampler number to use based on the current iteration
  unsigned int GetSamplerNumber(int iter);
  
  /// Abstract Samplers - for the MCMC SAEM
  std::vector<std::shared_ptr<AbstractSampler>> samplers_;

  /// Number of iterations to do for each sampler
  std::vector<unsigned int> number_of_iterations_per_sampler_;
  
  /// Candidate random variables
  std::shared_ptr<CandidateRandomVariables> candidate_rand_var_;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Abstract Model
  std::shared_ptr<AbstractModel> model_;

  /// Awesome Rea
  std::shared_ptr<Realizations> realizations_;

  /// Stochastic Sufficient Statistics used in the stochastic approximation step
  SufficientStatisticsVector stochastic_sufficient_stats_;

  /// Total number of iterations
  unsigned int max_iter_num_;

  /// Number of burn-in iterations
  unsigned int burnin_iter_num_;
  
  Algorithm(const Algorithm&);
  Algorithm& operator=(const Algorithm&);
};
