#ifndef _BlockedGibbsSampler_h
#define _BlockedGibbsSampler_h


#include "AbstractSampler.h"

class BlockedGibbsSampler : public AbstractSampler {

public:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // For each block; vector of Name + SubjectNumber
    typedef std::vector< std::tuple< std::string, unsigned int>> Block;
    

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
    virtual void InitializeSampler(const std::shared_ptr<MultiRealizations>& R);
    
    // Sample a new variable thanks to the sampler
    // The model cannot be constant because we modify some of its parameters (m_Orthonormal Basis for instance)
    virtual MultiRealizations Sample(const std::shared_ptr<MultiRealizations>& R, std::shared_ptr<AbstractModel>& M,
                                     const std::shared_ptr<Data>& D, int IterationNumber);


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Sample one block
    MultiRealizations OneBlockSample(int BlockNumber, const std::shared_ptr<MultiRealizations>& R, 
                                     std::shared_ptr<AbstractModel>& M,
                                     const std::shared_ptr<Data>& D, int IterationNumber);
    
    /// Check if all the random variables are from one individual
    bool IndividualRandomVariables(Block B);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Candidates random variables, corresponding to those in the Model
    CandidateRandomVariables m_CandidateRandomVariables;
    
    /// Blocks of the sampler
    std::vector<Block> m_Blocks;
    
    /// Last likelihood computed : reused 
    double m_LastLikelihoodComputed = 0;
    
    /// Names of the random variables that are individual
    std::vector<std::string> m_IndividualRandomVariables;
    
};


#endif //_BlockedGibbsSampler_h
