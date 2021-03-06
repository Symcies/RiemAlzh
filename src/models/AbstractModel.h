#ifndef _AbstractModel_h
#define _AbstractModel_h

#include <algorithm>
#include <iostream>
#include <fstream>
#include <functional>
#include <math.h>
#include <memory>
#include <unordered_map>


#include "Observations.h"
#include "ModelSettings.h"
#include "MultiRandomVariables.h"
#include "Realizations.h"
#include "ReadData.h"
#include "LinearAlgebra.h"
#include "TestAssert.h"
#include "UniformRandomVariable.h"
#include "AbstractRandomVariable.h"
#include "AbstractManifold.h"


class AbstractModel {
public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

      
  typedef typename LinearAlgebra<ScalarType>::MatrixType MatrixType;
  typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
  typedef typename std::unordered_map<std::string, int> StringIntHash;
  typedef std::vector<std::pair<std::string,  unsigned int>> MiniBlock;
  typedef std::pair<int, MiniBlock> SamplerBlock;
  typedef std::vector<VectorType> SufficientStatisticsVector;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  std::shared_ptr< AbstractRandomVariable > GetRandomVariable(std::string name) const;
  std::shared_ptr< AbstractRandomVariable > GetRandomVariable(int Key) const;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initialize the model
  virtual void Initialize(const Observations& Obs) = 0;
  
  /// Initialize the variance of the proposition distribution
  virtual ScalarType InitializePropositionDistributionVariance(std::string Name) const = 0;
      
  /// Update parameters ; some model-specifid private members need to be initilize, m_Orthogonal Basis for instance
  /// This update can depend on the parameter that has changed, provided by the Name argument
  virtual void UpdateModel(const Realizations& AR, int Type, const std::vector<std::string> Names = {"All"}) = 0;

  /// Update the sufficient statistics according to the model variables / parameters 
  virtual SufficientStatisticsVector GetSufficientStatistics(const Realizations& AR, const Observations& Obs) = 0;
  
  /// Update the fixed effects thanks to the approximation step of the algorithm
  virtual void UpdateRandomVariables(const SufficientStatisticsVector& StochSufficientStatistics) = 0;
  
  /// Compute the log likelihood of the model
  /// Using the log likelihood may have computational reason - for instance when the likelihood is too small
  virtual double ComputeLogLikelihood(const Observations &Obs)= 0;
  
  /// Compute the log likelihood of the model for a particular individual
  virtual double ComputeIndividualLogLikelihood(const IndividualObservations& Obs ,const int SubjectNumber) = 0;
  
  /// Simulate data according to the model
  virtual Observations SimulateData(io::DataSettings& DS) = 0;

  /// Simulate some random variable realizations
  Realizations SimulateRealizations();
  
  /// Define the sampler block used in the gibbs sampler (should it be here?)
  virtual std::vector<SamplerBlock> GetSamplerBlocks() const = 0;

  /// PROBABLY TO ERASE
  /// Compute the parallel curve
  virtual VectorType ComputeParallelCurve(int SubjectNumber, int ObservationNumber) = 0;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Outputs
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /// Compute Outputs
  virtual void DisplayOutputs(const Realizations& AR) = 0;
  
  /// Save the data into a file
  virtual void SaveData(unsigned int IterationNumber, const Realizations& R) = 0;
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initialize the true parameters to simulate data according to it - these parameters are unknown to the algo
  virtual void InitializeFakeRandomVariables() = 0;


protected:
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Riemanian manifold
  std::shared_ptr< AbstractManifold > m_Manifold;

  /// Random variables 
  MultiRandomVariables m_RandomVariables;
  
  /// Number of realizations per random variables used in the model
  StringIntHash m_RealizationsPerRandomVariable;
  
  /// Output file
  std::ofstream m_OutputParameters;
  
  /// Sum of the observations - corresponds to the first sufficient statistic
  double m_SumObservations;
  
  /// Total number of observations throughout all the individuals
  double m_NbTotalOfObservations;
  
  /// Total number of subjects
  unsigned int m_NumberOfSubjects;
  
  /// Dimension of the manifold
  double m_ManifoldDimension;

};


#endif //_AbstractModel_h
