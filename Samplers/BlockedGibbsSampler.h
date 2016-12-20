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

    BlockedGibbsSampler();


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

   
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the sampler
    virtual void InitializeSampler(const std::shared_ptr<Realizations>& R);
    
    // Sample a new variable thanks to the sampler
    // The model cannot be constant because we modify some of its parameters (m_Orthonormal Basis for instance)
    virtual Realizations Sample(const std::shared_ptr<Realizations>& R, std::shared_ptr<AbstractModel>& M,
                                     const std::shared_ptr<Data>& D, int IterationNumber);


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Sample one block
    Realizations OneBlockSample(int BlockNumber, const std::shared_ptr<Realizations>& R, 
                                     std::shared_ptr<AbstractModel>& M,
                                     const std::shared_ptr<Data>& D, int IterationNumber);
    
    
    /// Check if all the random variables are from one individual
    int TypeRandomVariables(Block B);
    
    /// Compute likelihood based on the block type
    VectorType ComputeLogLikelihood(int Type, const std::shared_ptr<Realizations> R, 
                                    const std::shared_ptr<AbstractModel> M, const std::shared_ptr<Data> D);
    
    /// Get previously computed log likelihood
    double GetPreviousLogLikelihood(int Type, const std::shared_ptr<AbstractModel> M, 
                                    const std::shared_ptr<Realizations> R, const std::shared_ptr<Data> D);
    
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
    
    
};


#endif //_BlockedGibbsSampler_h
