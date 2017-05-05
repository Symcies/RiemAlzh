#pragma once

#include <utility>

#include "AbstractModel.h"
#include "GaussianModel.h"

class GaussianMixtureModel: public AbstractModel {
public:
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  typedef std::unordered_map<int, int> IntClassHash;
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  GaussianMixtureModel(io::ModelSettings& model_settings);
  ~GaussianMixtureModel();

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  virtual std::shared_ptr< AbstractRandomVariable > GetRandomVariable(int key) const;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initialize the model
  virtual void Initialize(const Observations& obs);
  
  /// Initialize the model in case of validation data
  virtual void InitializeValidationDataParameters(const io::SimulatedDataSettings& data_settings, const io::ModelSettings& model_settings);
  
  /// Initialize the variance of the proposition distribution
  //virtual ScalarType InitializePropositionDistributionVariance(std::string name) const;
  
  /// Update the model parameters != random variables parameters
  virtual void UpdateModel(const Realizations& reals, const MiniBlock& block_info, const std::vector<std::string> names = {"All"});

  /// Update the sufficient statistics according to the model variables / parameters
  virtual SufficientStatisticsVector GetSufficientStatistics(const Realizations& reals, const Observations& obs);

  /// Update the fixed effects thanks to the approximation step of the algorithm
  virtual void UpdateRandomVariables(const SufficientStatisticsVector& stoch_sufficient_stats);

  /// Simulate data according to the model
  virtual Observations SimulateData(io::SimulatedDataSettings& data_settings);

  /// Define the sampler block used in the gibbs sampler (should it be here?)
  virtual std::vector<MiniBlock> GetSamplerBlocks() const;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Log-likelihood related method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /// Initialize the loglikelihood vector of the model
  virtual void InitializeLogLikelihood(const Observations& obs);
  
  /// Compute the log likelihood of the model
  virtual VectorType ComputeLogLikelihood(const Observations &obs, const MiniBlock& block_info);

  /// Compute the log likelihood of the model for a particular individual
  virtual ScalarType ComputeIndividualLogLikelihood(const IndividualObservations& obs ,const int subject_num);
  
  /// Get the previous loglikelihood computed
  virtual ScalarType GetPreviousLogLikelihood(const MiniBlock& block_info);
  
  /// Update the previous loglikelihood computed
  virtual void SetPreviousLogLikelihood(VectorType& log_likelihood, const MiniBlock& block_info);
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Outputs
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Compute Outputs
  virtual void DisplayOutputs(const Realizations& reals);

  /// Save the data into a file
  virtual void SaveCurrentState(unsigned int IterationNumber, const Realizations& reals);
  
  /// Save the final parameters and realizations into a file
  virtual void SaveFinalState(const Realizations& reals);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Overwriten method(s):
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  virtual Realizations SimulateRealizations();
  
  
 private:
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /// Based on the block_info, get the class of the block --> There should be only one class per block
  int GetClassNumber(const MiniBlock& block_info);
  
  /// Get the realizations of a specific class
  Realizations GetClassRealizations(const unsigned int class_number, const Realizations& reals);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Number of classes / sub-populations
  unsigned int number_of_classes_;
  
  /// Pointers to all the classes
  std::vector<GaussianModel> models_;
  
  /// Probabilities for each class
  VectorType class_probabilities;
  
  /// vector of sufficient statistics size
  std::vector<unsigned int> suff_stat_sizes_;
  
  /// Matrix of log-likelihoods : row -> class ; column -> individual
  MatrixType last_loglikelihood_;
  
  IntClassHash    key_to_class_;
};

