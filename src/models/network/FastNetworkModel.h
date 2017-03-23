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

  /// Initialize the variance of the proposition distribution
  virtual ScalarType InitializePropositionDistributionVariance(std::string Name) const;

  /// Update the model parameters != random variables parameters
  virtual void UpdateModel(const Realizations& R, int Type, const std::vector<std::string> Names = {"All"});

  /// Get the sufficient statistics of the model
  virtual SufficientStatisticsVector GetSufficientStatistics(const Realizations& R, const Observations& Obs);

  /// Update the random variables <=> the parameters of the model
  virtual void UpdateRandomVariables(const SufficientStatisticsVector& SS);

  /// Compute the log likelihood of the model
  virtual double ComputeLogLikelihood(const Observations &Obs);

  /// Compute the log likelihood of the model for a given subject
  virtual double ComputeIndividualLogLikelihood(const IndividualObservations& Obs, const int SubjectNumber);

  /// Simulate data according to the model and the parameters
  virtual Observations SimulateData(io::DataSettings& DS);

  /// Define the sampler block used in the gibbs sampler (should it be here?)
  virtual std::vector<SamplerBlock> GetSamplerBlocks() const ;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Outputs
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Compute outputs
  virtual void DisplayOutputs(const Realizations& R);

  /// Save the data into a file
  virtual void SaveData(unsigned int IterationNumber, const Realizations& R);

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initialize the true parameters to simulate data according to it - these parameters are unknown to the algo
  virtual void InitializeFakeRandomVariables();

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
  std::vector<VectorType> indiv_obs_date_;

  /// Time reparametrization of each individual
  std::vector<VectorType> indiv_time_points_;

  /// Block1 corresponds to p0 * exp(Delta)
  VectorType block1_;

  /// Block2 corresponds to vu_k / p0
  VectorType block2_;


private:
  FastNetworkModel(const FastNetworkModel &);
  FastNetworkModel& operator=(const FastNetworkModel &);
};
