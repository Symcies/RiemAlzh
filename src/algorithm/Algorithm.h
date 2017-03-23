#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#include "AbstractModel.h"
#include "AbstractSampler.h"
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
  // Constructor(stat_vector) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  Algorithm(io::AlgorithmSettings& settings);
  ~Algorithm();

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Encapsulation method(stat_vector) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  inline void SetModel(const std::shared_ptr<AbstractModel>& m ) { model_ = m; }

  inline void SetSampler(std::shared_ptr<AbstractSampler> stat_vector) { sampler_ = stat_vector; }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Other method(stat_vector) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compute the MCMC SAEM algorithm
    void ComputeMCMCSAEM(const Observations& obs);
    /// Compute the MCMC SAEM algorithm
    void IterationMCMCSAEM(const Observations& obs, int iter);

private:
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Method(stat_vector) of the MCMC SAEM:
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initialize the stochastic approximation
  void InitializeStochasticSufficientStatistics(const Observations& obs);

  /// Initialize sampler
  void InitializeSampler();

  /// Initialize Manifold
  void InitializeModel(const Observations& obs);

  /// Compute the simulation step : Gibbs Sampling
  void ComputeSimulationStep(const Observations& obs, int iter);

  /// Compute the stochastic coefficient
  void ComputeStochasticApproximation(SufficientStatisticsVector& stat_vector, int iter);

  /// Compute the decreasing step size of the approximation step
  double DecreasingStepSize(int iter);

  bool IsOutputIteration(int iter);
  bool IsDataSaveIteration(int iter);

  void DisplayIterations(int iter);
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Output(stat_vector)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Number of iterations to wait till next output display
  //TODO: should be const?
  unsigned int output_iter_;

  /// Number of iterations to wait till next data saving
  //TODO: should be const?
  unsigned int data_save_iter_;

  /// Acceptance ratios
  IntVectorHash acceptance_ratio_;

  /// Compute the acceptance ratio for each random variable
  void ComputeAcceptanceRatio(Realizations& prev_real_, int iter);

  /// Display acceptance ratio
  void DisplayAcceptanceRatio();

  /// Display Outputs
  void DisplayOutputs();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Attribute(stat_vector)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Abstract Model
  std::shared_ptr<AbstractModel> model_;

  /// Awesome Rea
  std::shared_ptr<Realizations> realizations_;

  /// Abstract Sampler - for the MCMC SAEM
  std::shared_ptr<AbstractSampler> sampler_;

  /// Stochastic Sufficient Statistics used in the stochastic approximation step
  SufficientStatisticsVector stochastic_sufficient_stats_;

  /// Total number of iterations
  unsigned int max_iter_num_;

  /// Number of burn-in iterations
  unsigned int burnin_iter_num_;

  /// Number of iterations done by the MCMC-SAEM
  unsigned int iter_count_ = 0;

  Algorithm(const Algorithm&);
  Algorithm& operator=(const Algorithm&);
};
