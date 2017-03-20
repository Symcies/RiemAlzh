#pragma once

#include "AbstractModel.h"

class UnivariateModel : public AbstractModel {
 public:
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  UnivariateModel(io::ModelSettings& model_settings);
  ~UnivariateModel();

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initialize the model
  virtual void Initialize(const Observations& obs);

  /// Initialize the variance of the proposition distribution
  virtual ScalarType InitializePropositionDistributionVariance(std::string name) const;

  /// Update the model parameters != random variables parameters
  virtual void UpdateModel(const Realizations& reals, int type, const std::vector<std::string> names = {"All"});

  /// Update the sufficient statistics according to the model variables / parameters
  virtual SufficientStatisticsVector GetSufficientStatistics(const Realizations& reals, const Observations& obs);

  /// Update the fixed effects thanks to the approximation step of the algorithm
  virtual void UpdateRandomVariables(const SufficientStatisticsVector& stoch_sufficient_stats);

  /// Compute the log likelihood of the model
  virtual ScalarType ComputeLogLikelihood(const Observations &obs);

  /// Compute the log likelihood of the model for a particular individual
  virtual ScalarType ComputeIndividualLogLikelihood(const IndividualObservations& obs, const int subjects_tot_num_);

  /// Simulate data according to the model
  virtual Observations SimulateData(io::DataSettings& data_settings);

  /// Define the sampler block used in the gibbs sampler (should it be here?)
  virtual std::vector<SamplerBlock> GetSamplerBlocks() const;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Outputs
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Compute Outputs
  virtual void DisplayOutputs(const Realizations& reals);

  /// Save the data into a file
  virtual void SaveData(unsigned int IterationNumber, const Realizations& reals);

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initialize the true parameters to simulate data according to it - these parameters are unknown to the algo
  virtual void InitializeFakeRandomVariables();

  /// Compute the parallel curve
  VectorType ComputeParallelCurve(int subjects_tot_num_, int ObservationNumber);


protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Compute the subject time points
  void ComputeSubjectTimePoint(const Realizations& reals, const int subjects_tot_num_ = -1);

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Noise model
  std::shared_ptr< GaussianRandomVariable > noise_;

  /// Real time of observation of each individual
  std::vector<VectorType> individual_obs_date_;

  /// Time reparametrization of each individual
  std::vector<VectorType> subj_time_points_;

  /// Attribute encoding for the position P
  ScalarType position_;

};
