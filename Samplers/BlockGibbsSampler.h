#ifndef _BlockGibbsSampler_h
#define _BlockGibbsSampler_h


#include "AbstractSampler.h"

class BlockGibbsSampler : public AbstractSampler {
public:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    typedef std::vector< std::pair< std::string, int>> Block;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    BlockGibbsSampler();
    ~BlockGibbsSampler();


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

   
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the sampler
    virtual void InitializeSampler(const std::shared_ptr<MultiRealizations>& R);
    
    // Sample a new variable thanks to the sampler
    virtual MultiRealizations Sample(const std::shared_ptr<MultiRealizations>& R, std::shared_ptr<AbstractModel>& M,
                                     const std::shared_ptr<Data>& D, int IterationNumber);
    
protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Sample one block 
    MultiRealizations BlockSample(int BlockNumber, const MultiRealizations& R, std::shared_ptr<AbstractModel>& M, 
                                  const std::shared_ptr<Data>& D, int IterationNumber);
    

    /// Compute likelihood
    VectorType ComputeLikelihood(const std::shared_ptr<MultiRealizations>& R, std::shared_ptr<AbstractModel>& M, 
                             const std::shared_ptr<Data>& D, int Type);
            
    /// Get the previously computed likelihood
    double GetPreviousLogLikelihood(int Type);
    
    
    /// Check the type of the random variables in the block --> In order to compute the loglikelihood
    int TypeRandomVariables(Block B);
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
    /// Blocks of the sampler
    std::vector<Block> m_Blocks;
    
    /// Last loglikelihood computed : each coordinate is an individual likelihood
    VectorType m_LogLikelihood;
    
    /// Candidate random variables
    CandidateRandomVariables m_CandidateRandomVariables;
    
    
};


#endif //_BlockGibbsSampler_h
