#pragma once

#include <cmath>
#include <iostream>
#include <fstream>
#include <memory>
#include <unordered_map>
#include <tuple>

#include "AbstractManifold.h"
#include "AbstractRandomVariable.h"
#include "SimulatedDataSettings.h"
#include "LinearAlgebra.h"
#include "ModelSettings.h"
#include "MultiRandomVariables.h"
#include "Observations.h"
#include "Realizations.h"
#include "ReadData.h"
#include "UniformRandomVariable.h"
#include "global.h"


class AbstractModel {
public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  typedef typename LinearAlgebra<ScalarType>::MatrixType MatrixType;
  typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
  typedef typename std::unordered_map<std::string, int> StringIntHash;
  typedef std::vector<std::tuple<int, std::string, int>> MiniBlock;
  typedef std::vector<VectorType> SufficientStatisticsVector;
  typedef std::unordered_map<std::string, std::pair<std::vector<double>, ScalarType>> InitialRVParameters;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  AbstractModel(){};
  virtual ~AbstractModel(){};

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  virtual std::shared_ptr< AbstractRandomVariable > GetRandomVariable(std::string name) const;
  virtual std::shared_ptr< AbstractRandomVariable > GetRandomVariable(int key) const;

  inline std::shared_ptr< AbstractManifold > GetManifold(){return manifold_;};
  inline MultiRandomVariables& GetRandomVariable(){return rand_var_;};
  inline StringIntHash GetAssoTable(){return asso_num_real_per_rand_var_;};
  inline double GetSumOfObservations(){return sum_obs_;};
  inline double GetNumberOfObservations(){return obs_tot_num_;};
  inline int GetNumberOfSubjects(){return subjects_tot_num_;};
  inline double GetManifoldDimension(){return manifold_dim_;};
  inline std::vector<std::string> GetAcceptanceRatioToDisplay() { return acceptance_ratio_to_display_; }
  
  inline void SetRandomVariableParameters(const InitialRVParameters& params) { rv_params_ = params; } 

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initialize the model
  virtual void Initialize(const Observations& obs) = 0;
  
  /// Initialize the model in case of validation data
  virtual void InitializeValidationDataParameters(const io::SimulatedDataSettings& data_settings, const io::ModelSettings& model_settings) = 0;
  
  /// Initialize the variance of the proposition distribution
  virtual ScalarType InitializePropositionDistributionVariance(std::string name) const;

  /// Update parameters ; some model-specifid private members need to be initilize, m_Orthogonal Basis for instance
  /// This update can depend on the parameter that has changed, provided by the name argument
  virtual void UpdateModel(const Realizations& real, const MiniBlock& block_info, const std::vector<std::string> names = {"All"}) = 0;

  /// Update the sufficient statistics according to the model variables / parameters
  virtual SufficientStatisticsVector GetSufficientStatistics(const Realizations& real, const Observations& obs) = 0;

  /// Update the fixed effects thanks to the approximation step of the algorithm
  virtual void UpdateRandomVariables(const SufficientStatisticsVector& stoch_sufficient_stats) = 0;
  
  /// Simulate data according to the model
  virtual Observations SimulateData(io::SimulatedDataSettings& data_settings) = 0;

  /// Simulate some random variable realizations
  virtual Realizations SimulateRealizations();

  /// Define the sampler block used in the gibbs sampler (should it be here?)
  virtual std::vector<MiniBlock> GetSamplerBlocks() const = 0;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Log-likelihood related method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /// Initialize the loglikelihood vector of the model
  virtual void InitializeLogLikelihood(const Observations& obs) = 0;
  
  /// Compute the log likelihood of the model
  virtual VectorType ComputeLogLikelihood(const Observations &obs, const MiniBlock& block_info)= 0;

  /// Compute the log likelihood of the model for a particular individual
  virtual ScalarType ComputeIndividualLogLikelihood(const IndividualObservations& obs, const int subjects_tot_num_) = 0;
  
  /// Get the previous loglikelihood computed
  virtual ScalarType GetPreviousLogLikelihood(const MiniBlock& block_info) = 0;
  
  /// Update the previous loglikelihood computed
  virtual void SetPreviousLogLikelihood(VectorType& log_likelihood, const MiniBlock& block_info) = 0;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Outputs
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Compute Outputs
  virtual void DisplayOutputs(const Realizations& reals) = 0;

  /// Save the current parameters and realizations into a file
  virtual void SaveCurrentState(unsigned int iter_num, const Realizations& reals) = 0;
  
  /// Save the final parameters and realizations into a file
  virtual void SaveFinalState(const Realizations& reals, const Observations& obs) = 0;
  virtual void SavePopulationFile() = 0;
  virtual void SaveIndividualsFile(const Realizations &reals, const Observations &obs) = 0;

protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Riemanian manifold
  std::shared_ptr< AbstractManifold > manifold_;

  /// Random variables
  MultiRandomVariables rand_var_;
  
  /// Number of realizations per random variables used in the model
  StringIntHash asso_num_real_per_rand_var_;

  /// Random variable proposition distribution
  std::unordered_map<std::string, ScalarType> proposition_distribution_variance_;

  /// Initial Random variable parameters : type of RV, mean, variance, proposition distribution, ...
  InitialRVParameters rv_params_;
  
  /// Output file
  std::string output_file_name_;

  /// Sum of the observations - corresponds to the first sufficient statistic
  double sum_obs_;

  /// Total number of observations throughout all the individuals
  double obs_tot_num_;

  /// Total number of subjects
  unsigned int subjects_tot_num_;

  /// Dimension of the manifold
  double manifold_dim_;
  
  /// Name of the random variable acceptance ratio to displau
  std::vector<std::string> acceptance_ratio_to_display_;

private:
  AbstractModel(const AbstractModel &);
  AbstractModel& operator=(const AbstractModel &);

};
