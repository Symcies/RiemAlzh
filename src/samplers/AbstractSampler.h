#ifndef _AbstractSampler_h
#define _AbstractSampler_h


#include <memory>
#include <algorithm>
#include <string>
#include <map>

#include "GaussianRandomVariable.h"
#include "AbstractRandomVariable.h"
#include "AbstractModel.h"
#include "CandidateRandomVariables.h"
#include "Realizations.h"

class AbstractSampler {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
        
    // For each block; vector of Name + SubjectNumber
    typedef std::vector< std::tuple< std::string, int>> Block;
    
    
    typedef typename LinearAlgebra<ScalarType>::MatrixType MatrixType;
    typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
    
    typedef std::unordered_map<std::string, unsigned int> MiniBlock;
    typedef std::pair<int, MiniBlock> SamplerBlock;
    
    typedef std::vector< std::vector< std::pair< VectorType, double> > > Data;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    AbstractSampler();
    ~AbstractSampler();


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

   
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the sampler
    virtual void InitializeSampler(Realizations& R, AbstractModel &M, const Data& D) = 0;
    
    /// Sample new realizations of the model random variables
    virtual void Sample(Realizations& R, AbstractModel& M, const Data& D) = 0;


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compute the decreasing step size of the adaptive variances
    double DecreasingStepSize(int Iteration, int NoMemoryTime);
    
    /// Update the variance of the gaussian proposition distribution
    void UpdatePropositionDistributionVariance(GaussianRandomVariable& GRV, double Ratio, int Iteration);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////
        
    /// Blocks of the sampler
    std::vector<SamplerBlock> m_Blocks;
    
    /// Sampling time without memory
    unsigned int m_MemorylessSamplingTime = 10000;
    
    /// Acceptation ratio 
    double m_ExpectedAcceptanceRatio = 0.301;
    
    /// Candidates random variables, corresponding to those in the Model
    CandidateRandomVariables m_CandidateRandomVariables;
    

};


#endif //_AbstractSampler_h
