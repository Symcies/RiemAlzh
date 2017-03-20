#ifndef _AbstractModel_h
#define _AbstractModel_h

#include <cmath>
#include <iostream>
#include <fstream>
#include <memory>
#include <unordered_map>

#include "AbstractManifold.h"
#include "AbstractRandomVariable.h"
#include "DataSettings.h"
#include "LinearAlgebra.h"
#include "ModelSettings.h"
#include "MultiRandomVariables.h"
#include "Observations.h"
#include "Realizations.h"
#include "ReadData.h"
#include "UniformRandomVariable.h"


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
  std::shared_ptr< AbstractRandomVariable > GetRandomVariable(int key) const;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initialize the model
  virtual void Initialize(const Observations& obs) = 0;

  /// Initialize the variance of the proposition distribution
  virtual ScalarType InitializePropositionDistributionVariance(std::string name) const = 0;

  /// Update parameters ; some model-specifid private members need to be initilize, m_Orthogonal Basis for instance
  /// This update can depend on the parameter that has changed, provided by the name argument
  virtual void UpdateModel(const Realizations& real, int type, const std::vector<std::string> names = {"All"}) = 0;

  /// Update the sufficient statistics according to the model variables / parameters
  virtual SufficientStatisticsVector GetSufficientStatistics(const Realizations& real, const Observations& obs) = 0;

  /// Update the fixed effects thanks to the approximation step of the algorithm
  virtual void UpdateRandomVariables(const SufficientStatisticsVector& stoch_sufficient_stats) = 0;

  /// Compute the log likelihood of the model
  /// Using the log likelihood may have computational reason - for instance when the likelihood is too small
  virtual double ComputeLogLikelihood(const Observations &obs)= 0;

  /// Compute the log likelihood of the model for a particular individual
  virtual double ComputeIndividualLogLikelihood(const IndividualObservations& obs ,const int subjects_tot_num_) = 0;

  /// Simulate data according to the model
  virtual Observations SimulateData(io::DataSettings& data_settings) = 0;

  /// Simulate some random variable realizations
  Realizations SimulateRealizations();

  /// Define the sampler block used in the gibbs sampler (should it be here?)
  virtual std::vector<SamplerBlock> GetSamplerBlocks() const = 0;

  /// PROBABLY TO ERASE
  /// Compute the parallel curve
  virtual VectorType ComputeParallelCurve(int subjects_tot_num_, int obs_num) = 0;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Outputs
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Compute Outputs
  virtual void DisplayOutputs(const Realizations& reals) = 0;

  /// Save the data into a file
  virtual void SaveData(unsigned int iter_num, const Realizations& reals) = 0;


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
  std::shared_ptr< AbstractManifold > manifold_;

  /// Random variables
  MultiRandomVariables rand_var_;

  /// Number of realizations per random variables used in the model
  //TODO: why StringIntHash
  StringIntHash asso_num_real_per_rand_var_;

  /// Output file
  std::ofstream out_params_;

  /// Sum of the observations - corresponds to the first sufficient statistic
  double sum_obs_;

  /// Total number of observations throughout all the individuals
  double obs_tot_num_;

  /// Total number of subjects
  unsigned int subjects_tot_num_;

  /// Dimension of the manifold
  double manifold_dim_;

};


#endif //_AbstractModel_h
