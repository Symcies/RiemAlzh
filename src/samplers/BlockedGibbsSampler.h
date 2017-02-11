#ifndef _BlockedGibbsSampler_h
#define _BlockedGibbsSampler_h

#include "AbstractSampler.h"
#include <algorithm>

class BlockedGibbsSampler : public AbstractSampler {

public:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // For each block; vector of Name + SubjectNumber
    typedef std::vector< std::tuple< std::string, int>> Block;
    

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
    virtual void InitializeSampler(const Realizations &R, AbstractModel &M, const Data& D);
    
    /// Sample new realizations
    virtual Realizations Sample(Realizations& R, AbstractModel& M, const Data& D, int IterationNumber);


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Sample one block
    void OneBlockSample(int BlockNumber, Realizations& R, AbstractModel& M, const Data& D, int IterationNumber);
    
        /// Compute Prior part of the ratio while updating the realization
    ScalarType ComputePriorRatioAndUpdateRealizations(Realizations& R, const AbstractModel& M, const Block& Variables);
    
    /// Check if all the random variables are from one individual
    int TypeRandomVariables(Block& B);
    
    /// Compute likelihood based on the block type
    VectorType ComputeLogLikelihood(int Type, const Realizations& R, AbstractModel& M, const Data& D);
    
    /// Get previously computed log likelihood
    double GetPreviousLogLikelihood(int Type, AbstractModel& M, const Realizations& R, const Data& D);
    
    /// Update the last log likelihood computed
    void UpdateLastLogLikelihood(int Type, VectorType& ComputedLogLikelihood);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Candidates random variables, corresponding to those in the Model
    CandidateRandomVariables m_CandidateRandomVariables;
    
    /// Blocks of the sampler
    std::vector<Block> m_Blocks;
    
    /// Last likelihood computed : reused 
    VectorType m_LastLikelihoodComputed;
    
        /// Parameters updated in the block
    std::vector<std::string> m_CurrentBlockParameters;
    std::vector<int> m_CurrentBlockParametersBIS;
    
    /// Parameters to recover back in the realization
    std::unordered_map<std::string, std::pair<unsigned int, ScalarType >> m_RecoverParameters;
    
    
};


#endif //_BlockedGibbsSampler_h
