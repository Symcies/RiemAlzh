#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <math.h>

#include "Observations.h"
#include "AlgorithmSettings.h"
#include "AbstractModel.h"
#include "AbstractSampler.h"
#include "CandidateRandomVariables.h"


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

  Algorithm(io::AlgorithmSettings& Settings);
  

  ~Algorithm();

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  inline void SetModel(const std::shared_ptr<AbstractModel>& M ) { m_Model = M; }

  inline void SetSampler(std::shared_ptr<AbstractSampler> S) { m_Sampler = S; }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Compute the MCMC SAEM algorithm
  void ComputeMCMCSAEM(const Observations& Obs);


protected:
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Method(s) of the MCMC SAEM:
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initialize the stochastic approximation
  void InitializeStochasticSufficientStatistics(const Observations& Obs);
      
  /// Initialize sampler
  void InitializeSampler();
  
  /// Initialize Manifold
  void InitializeModel(const Observations& Obs);

  /// Compute the simulation step : Gibbs Sampling
  void ComputeSimulationStep(const Observations& Obs);

  /// Compute the stochastic coefficient 
  void ComputeStochasticApproximation(SufficientStatisticsVector& S);

  /// Compute the decreasing step size of the approximation step
  double DecreasingStepSize();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Output(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Number of iterations to wait till next output display
  unsigned int m_CounterToDisplayOutputs;
  
  /// Number of iterations to wait till next data saving
  unsigned int m_CounterToSaveData;
  
  /// Acceptance ratios
  IntVectorHash m_AcceptanceRatios;
  
  /// Compute the acceptance ratio for each random variable
  void ComputeAcceptanceRatio(Realizations& PreviousReals);
  
  /// Display acceptance ratio
  void DisplayAcceptanceRatio();
  
  /// Display Outputs
  void DisplayOutputs();
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Abstract Model
  std::shared_ptr<AbstractModel> m_Model;
      
  /// Awesome Rea
  std::shared_ptr<Realizations> m_Realizations;

  /// Abstract Sampler - for the MCMC SAEM 
  std::shared_ptr<AbstractSampler> m_Sampler;

  /// Stochastic Sufficient Statistics used in the stochastic approximation step
  SufficientStatisticsVector m_StochasticSufficientStatistics;
  
  /// Total number of iterations
  unsigned int m_MaxNumberOfIterations;
  
  /// Number of burn-in iterations
  unsigned int m_BurnIn;
  
  /// Number of iterations done by the MCMC-SAEM
  unsigned int m_IterationCounter = 0;
};

