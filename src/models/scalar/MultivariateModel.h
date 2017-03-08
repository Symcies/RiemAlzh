#pragma once

#include "AbstractModel.h"

class MultivariateModel : public AbstractModel {
public:
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  MultivariateModel(io::ModelSettings& MS);
  ~MultivariateModel();
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initialize the model
  virtual void Initialize(const Observations& Obs);
  
  /// Initialize the variance of the proposition distribution
  virtual ScalarType InitializePropositionDistributionVariance(std::string Name) const;
  
  /// Update the model parameters != random variables parameters
  virtual void UpdateModel(const Realizations& R, int Type, const std::vector<std::string> Names = {"All"});

  /// Update the sufficient statistics according to the model variables / parameters 
  virtual SufficientStatisticsVector GetSufficientStatistics(const Realizations& R, const Observations& Obs);
  
  /// Update the fixed effects thanks to the approximation step of the algorithm
  virtual void UpdateRandomVariables(const SufficientStatisticsVector& StochSufficientStatistics);
  
  /// Compute the log likelihood of the model
  virtual double ComputeLogLikelihood(const Observations &Obs);
  
  /// Compute the log likelihood of the model for a particular individual
  virtual double ComputeIndividualLogLikelihood(const IndividualObservations& Obs, const int SubjectNumber);

  /// Simulate data according to the model
  virtual Observations SimulateData(io::DataSettings& DS);
  
  /// Define the sampler block used in the gibbs sampler (should it be here?)
  virtual std::vector<SamplerBlock> GetSamplerBlocks() const;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Outputs
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /// Compute Outputs
  virtual void DisplayOutputs(const Realizations& R);
  
  /// Save the data into a file
  virtual void SaveData(unsigned int IterationNumber, const Realizations& R);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initialize the true parameters to simulate data according to it - these parameters are unknown to the algo
  virtual void InitializeFakeRandomVariables();
  
  /// Probably to erase
  /// Compute the parallel curve
  VectorType ComputeParallelCurve(int SubjectNumber, int ObservationNumber);
  
  
protected:
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Compute the subjects time points
  void ComputeSubjectTimePoint(const Realizations& R, const int SubjectNumber = -1);
  
  /// Compute the delta
  void ComputeDeltas(const Realizations& R);
  
  /// Compute Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns>
  void ComputeOrthonormalBasis(); 

  /// Compute the A Matrix used to get the space shifts
  void ComputeAMatrix(const Realizations& R);

  /// Compute the space shifts
  void ComputeSpaceShifts(const Realizations& R);
  
  /// Compute block 1 (1/p0 -1)
  void ComputeBlock(const Realizations& R);
  

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Noise model
  std::shared_ptr< GaussianRandomVariable > m_Noise;
  
  /// Number of independent components
  unsigned int m_NbIndependentSources;
  
  /// Realization G encoding for the position p0 -> G = 1/P0 - 1
  ScalarType m_G;
  
  /// Realization encoding for the temporal shifts delta
  VectorType m_Deltas;
    
  /// Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns> (Basis orthogonal to gamma0_deriv(T0)
  MatrixType m_OrthogonalBasis;
  
  /// A Matrix vec<A1, ..., A(N)> where Ai is vec<Ns> (Ai is a column)
  MatrixType m_AMatrix;

  /// Space shifts w(i) of the model
  MatrixType m_SpaceShifts;
  
  /// Real time of observation of each individual
  std::vector<VectorType> m_IndividualObservationDate;
  
  /// Time reparametrization of each individual
  std::vector<VectorType> m_SubjectTimePoints;
      
  /// Block1 corresponds to p0 * exp(Delta)
  VectorType m_Block;
  
};
