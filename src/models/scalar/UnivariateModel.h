#ifndef _UnivariateModel_h
#define _UnivariateModel_h

#include "AbstractModel.h"

class UnivariateModel : public AbstractModel {
 public:
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  UnivariateModel(io::ModelSettings& MS);
  ~UnivariateModel();
  
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
  virtual ScalarType ComputeLogLikelihood(const Observations &Obs);
  
  /// Compute the log likelihood of the model for a particular individual
  virtual ScalarType ComputeIndividualLogLikelihood(const IndividualObservations& Obs, const int SubjectNumber);

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
  
  /// Compute the parallel curve
  VectorType ComputeParallelCurve(int SubjectNumber, int ObservationNumber);
  

protected:
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Compute the subject time points
  void ComputeSubjectTimePoint(const Realizations& R, const int SubjectNumber = -1);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Noise model
  std::shared_ptr< GaussianRandomVariable > m_Noise;
  
  /// Real time of observation of each individual
  std::vector<VectorType> m_IndividualObservationDate;
  
  /// Time reparametrization of each individual
  std::vector<VectorType> m_SubjectTimePoints;
  
  /// Attribute encoding for the position P
  ScalarType m_P;
  
};

#endif //_UnivariateModel_h
