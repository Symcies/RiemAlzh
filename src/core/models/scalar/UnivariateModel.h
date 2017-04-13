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
  
  /// Initialize the model in case of validation data
  virtual void InitializeValidationDataParameters(const io::SimulatedDataSettings& data_settings, const io::ModelSettings& model_settings);
  
  
  /// Update the model parameters != random variables parameters
  virtual void UpdateModel(const Realizations& reals, const MiniBlock& block_info, const std::vector<std::string> names = {"All"});

  /// Update the sufficient statistics according to the model variables / parameters
  virtual SufficientStatisticsVector GetSufficientStatistics(const Realizations& reals, const Observations& obs);

  /// Update the fixed effects thanks to the approximation step of the algorithm
  virtual void UpdateRandomVariables(const SufficientStatisticsVector& stoch_sufficient_stats);

  /// Simulate data according to the model
  virtual Observations SimulateData(io::SimulatedDataSettings& data_settings, bool need_init = false);

  /// Define the sampler block used in the gibbs sampler (should it be here?)
  virtual std::vector<MiniBlock> GetSamplerBlocks() const;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Log-likelihood related method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /// Compute the log likelihood of the model
  virtual VectorType ComputeLogLikelihood(const Observations &obs, const MiniBlock& block_info);

  /// Compute the log likelihood of the model for a particular individual
  virtual ScalarType ComputeIndividualLogLikelihood(const IndividualObservations& obs ,const int subjects_tot_num_);
  
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
  virtual void SaveData(unsigned int IterationNumber, const Realizations& reals);
  
  inline std::vector<VectorType> GetIndObsDate(){return individual_obs_date_;};
  inline std::vector<VectorType> GetSubjTimePoints(){return individual_time_points_;};

private:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  UnivariateModel(const UnivariateModel &);
  UnivariateModel& operator=(const UnivariateModel &);

  /// Compute the subject time points
  void ComputeSubjectTimePoint(const Realizations& reals, const int subjects_tot_num_ = -1);

  /// Compute the parallel curve
  VectorType ComputeParallelCurve(int subjects_tot_num_, int ObservationNumber);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Noise model
  std::shared_ptr< GaussianRandomVariable > noise_;

  /// Real time of observation of each individual
  std::vector<VectorType> individual_obs_date_;

  /// Time reparametrization of each individual
  std::vector<VectorType> individual_time_points_;

  /// Attribute encoding for the position P
  ScalarType position_;
  
  /// Last log-likelihood computed - vector of individual
  VectorType last_loglikelihood_;



};
