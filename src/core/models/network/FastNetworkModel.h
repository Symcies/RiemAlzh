#pragma once

#include "AbstractModel.h"

class FastNetworkModel : public AbstractModel {
public:
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  FastNetworkModel(io::ModelSettings& MS);
  ~FastNetworkModel();

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initialize the model : random variables, interpolation matrix, parameters
  virtual void Initialize(const Observations& Obs);
  
  /// Initialize the model in case of validation data
  virtual void InitializeValidationDataParameters(const io::SimulatedDataSettings& data_settings, const io::ModelSettings& model_settings);
  
  /// Update the model parameters != random variables parameters
  virtual void UpdateModel(const Realizations& reals, const MiniBlock& block_info, const std::vector<std::string> names = {"All"});

  /// Get the sufficient statistics of the model
  virtual SufficientStatisticsVector GetSufficientStatistics(const Realizations& R, const Observations& Obs);

  /// Update the random variables <=> the parameters of the model
  virtual void UpdateRandomVariables(const SufficientStatisticsVector& SS);

  /// Simulate data according to the model and the parameters
  virtual Observations SimulateData(io::SimulatedDataSettings& DS);

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
  virtual ScalarType ComputeIndividualLogLikelihood(const IndividualObservations& obs ,const int subjects_tot_num_);
  
  /// Get the previous loglikelihood computed
  virtual ScalarType GetPreviousLogLikelihood(const MiniBlock& block_info);
  
  /// Update the previous loglikelihood computed
  virtual void SetPreviousLogLikelihood(VectorType& log_likelihood, const MiniBlock& block_info);
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Outputs
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Compute outputs
  virtual void DisplayOutputs(const Realizations& R);

  /// Save the data into a file
  virtual void SaveCurrentState(unsigned int IterationNumber, const Realizations& R);

  /// Save the final parameters and realizations into a file
  virtual void SaveFinalState(const Realizations& reals);

protected:
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Compute the subjects time points
  void ComputeSubjectTimePoint(const Realizations& R, const int SubjectNumber = -1);

  /// Compute the interpolation coefficients delta
  void ComputeDeltas(const Realizations& R);

  /// Compute the interpolation coefficients beta
  void ComputeNus(const Realizations& R);

  /// Compute Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns>
  void ComputeOrthonormalBasis();

  /// Compute the A Matrix used to get the space shifts
  void ComputeAMatrix(const Realizations& R);

  /// Compute the space shifts
  void ComputeSpaceShifts(const Realizations& R);

  /// Compute the block p0 * exp(delta_k)
  void ComputeBlock1();

  /// Compute the block nu_k / p0
  void ComputeBlock2();
  
  /// Get the type of the sampler block
  int GetType(const MiniBlock& block_info);

  /// Compute the parallel curve
  VectorType ComputeParallelCurve(int SubjectNumber, int ObservationNumber);

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////


  /// Noise model
  std::shared_ptr< GaussianRandomVariable > noise_;

  /// Number of control points
  unsigned int control_points_nb_;

  /// Number of independent components
  unsigned int indep_components_nb_;

  /// Kernel Matrix K
  MatrixType invert_kernel_matrix_;

  /// Interpolation Matrix to calculate any interpolation
  MatrixType interpolation_matrix_;

  /// Initial position of the model P0 = exp(R.at("P0")(0))
  double init_pos_;

  /// Interpolation coefficients of delta
  VectorType deltas_;

  /// Interpolation coefficients of beta
  VectorType nus_;

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
  VectorType block1_;

  /// Block2 corresponds to vu_k / p0
  VectorType block2_;
  
  /// Last log-likelihood computed - vector of individual
  VectorType last_loglikelihood_;


private:
  FastNetworkModel(const FastNetworkModel &);
  FastNetworkModel& operator=(const FastNetworkModel &);
};
