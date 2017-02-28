#ifndef _BlockedGibbsSampler_h
#define _BlockedGibbsSampler_h

#include "AbstractSampler.h"
#include <algorithm>

class BlockedGibbsSampler : public AbstractSampler {

public:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Constructor
    BlockedGibbsSampler();
    
    /// Constructor 2
    BlockedGibbsSampler(unsigned int MemorylessSamplingTime, double ExpectedAcceptanceRatio);
    
    /// Destructor
    ~BlockedGibbsSampler();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

   
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the sampler
    virtual void InitializeSampler(Realizations& R, AbstractModel &M);
    
    /// Sample new realizations
    virtual void Sample(Realizations& R, AbstractModel& M, const OldData& D, const Observations& Obs);


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Sample one block
    void OneBlockSample(int BlockNumber, Realizations& R, AbstractModel& M, const OldData& D, const Observations& Obs);
    
    /// Compute Prior part of the ratio while updating the realization
    ScalarType ComputePriorRatioAndUpdateRealizations(Realizations& R, const AbstractModel& M, const MiniBlock& Variables);
    
    
    /// Compute likelihood based on the block type
    VectorType ComputeLogLikelihood(AbstractModel& M, const OldData& D, const Observations& Obs);
    
    /// Get previously computed log likelihood
    double GetPreviousLogLikelihood();
    
    /// Update the last log likelihood computed
    void UpdateLastLogLikelihood(VectorType& ComputedLogLikelihood);
    
        /// Update the random variable
    void UpdateBlockRandomVariable(double AcceptanceRatio, const MiniBlock& Variables);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Current Block Type
    int m_CurrentBlockType;
    
    /// Current iteration Number
    int m_CurrentIteration;
    

    
    /// Last likelihood computed : reused 
    VectorType m_LastLikelihoodComputed;
    
    /// Parameters updated in the block
    std::vector<std::string> m_CurrentBlockParameters;
    std::vector<int> m_CurrentBlockParametersBIS;
    
    /// Parameters to recover back in the realization
    std::vector<std::tuple<std::string, unsigned int, ScalarType >> m_RecoverParameters;
    
    /// Uniform distrubution
    std::uniform_real_distribution<double> m_UniformDistribution;
    
};


#endif //_BlockedGibbsSampler_h
