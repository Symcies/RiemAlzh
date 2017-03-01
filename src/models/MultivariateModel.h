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
        
    /// Update parameters ; some model-specifid private members need to be initilize, m_Orthogonal Basis for instance
    /// This update can depend on the parameter that has changed, provided by the Name argument
    virtual void UpdateModel(const Realizations& AR, int Type, const std::vector<std::string> Names = {"All"});

    /// Update the sufficient statistics according to the model variables / parameters 
    virtual SufficientStatisticsVector GetSufficientStatistics(const Realizations& AR, const Observations& Obs);
    
    /// Update the fixed effects thanks to the approximation step of the algorithm
    virtual void UpdateRandomVariables(const SufficientStatisticsVector& StochSufficientStatistics);
    
    
    /// Compute the log likelihood of the model
    /// Using the log likelihood may have computational reason - for instance when the likelihood is too small
    virtual double ComputeLogLikelihood(const Observations &Obs);
    
    /// Compute the log likelihood of the model for a particular individual
    virtual double ComputeIndividualLogLikelihood(const IndividualObservations& Obs, const int SubjectNumber);
    
    /// Simulate data according to the model
    virtual Observations SimulateData(io::DataSettings& DS);
    
    
    /// Define the sampler block used in the gibbs sampler (should it be here?)
    virtual std::vector<SamplerBlock> GetSamplerBlocks() const;

    /// PROBABLY TO ERASE
    /// Compute the parallel curve
    virtual VectorType ComputeParallelCurve(int SubjectNumber, int ObservationNumber);
    
    

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Outputs
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Compute Outputs
    virtual void DisplayOutputs(const Realizations& AR);
    
    /// Save the data into a file
    virtual void SaveData(unsigned int IterationNumber, const Realizations& R);
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the true parameters to simulate data according to it - these parameters are unknown to the algo
    virtual void InitializeFakeRandomVariables();
    
protected:
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Noise model
    std::shared_ptr< GaussianRandomVariable > m_Noise;
    
    /// Number of independent components
    unsigned int m_NbIndependentSources;
    
};
