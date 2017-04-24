#pragma once

#include "AbstractModel.h"

class MultivariateModel : public AbstractModel {
public:
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  MultivariateModel(io::ModelSettings& model_settings);
  ~MultivariateModel();

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
  virtual Observations SimulateData(io::SimulatedDataSettings& data_settings);

  /// Define the sampler block used in the gibbs sampler (should it be here?)
  virtual std::vector<MiniBlock> GetSamplerBlocks() const;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Log-likelihood related method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /// Compute the log likelihood of the model
  virtual VectorType ComputeLogLikelihood(const Observations &obs, const MiniBlock& block_info);

  /// Compute the log likelihood of the model for a particular individual
  virtual ScalarType ComputeIndividualLogLikelihood(const IndividualObservations& obs ,const int subject_num_);  
  
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



private:
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  MultivariateModel(const MultivariateModel &);
  MultivariateModel& operator=(const MultivariateModel &);

  /// Compute the subjects time points
  void ComputeSubjectTimePoint(const Realizations& reals, const int subjects_num = -1);

  /// Compute the delta
  void ComputeDeltas(const Realizations& reals);

  /// Compute Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns>
  void ComputeOrthonormalBasis();

  /// Compute the A Matrix used to get the space shifts
  void ComputeAMatrix(const Realizations& reals);

  /// Compute the space shifts
  void ComputeSpaceShifts(const Realizations& reals);

  /// Compute block 1 (1/p0 -1)
  void ComputeBlock(const Realizations& reals);
  
  /// Get the type of a sampler block
  int GetType(const MiniBlock& block_info);

  /// Compute the parallel curve
  VectorType ComputeParallelCurve(int subjects_num, int obs_num);


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Noise model
  std::shared_ptr< GaussianRandomVariable > noise_;

  /// Number of independent components
  unsigned int indep_sources_num_;

  /// Realization G encoding for the position p0 -> G = 1/P0 - 1
  ScalarType g_;

  /// Realization encoding for the temporal shifts delta
  VectorType deltas_;

  /// Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns> (Basis orthogonal to gamma0_deriv(T0)
  MatrixType orthog_basis_;

  /// A Matrix vec<A1, ..., A(N)> where Ai is vec<Ns> (Ai is a column)
  MatrixType a_matrix_;

  /// Space shifts w(i) of the model
  MatrixType space_shifts_;

  /// Real time of observation of each individual
  std::vector<VectorType> individual_obs_date_;

  /// Time reparametrization of each individual
  std::vector<VectorType> individual_time_points_;

  /// Block1 corresponds to p0 * exp(Delta)
  VectorType block_;
  
  /// Last log-likelihood computed - vector of individual
  VectorType last_loglikelihood_;

};
