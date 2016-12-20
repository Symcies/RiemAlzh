#ifndef _AbstractSampler_h
#define _AbstractSampler_h


#include <memory>
#include <algorithm>
#include <string>
#include "../RandomVariables/GaussianRandomVariable.h"
#include "../RandomVariables/AbstractRandomVariable.h"
#include <map>
#include "../Models/AbstractModel.h"
#include "../Parameters/CandidateRandomVariables.h"

class AbstractSampler {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    typedef typename LinearAlgebra<ScalarType>::MatrixType MatrixType;
    typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
    
    typedef std::vector< std::vector< std::pair< VectorType, double> > > Data;
    typedef std::map<std::string, VectorType> Realizations;

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
    virtual void InitializeSampler(const std::shared_ptr<Realizations>& R) = 0;
    
    // Sample a new variable thanks to the sampler
    // The model cannot be constant because we modify some of its parameters (m_Orthonormal Basis for instance)
    virtual Realizations Sample(const std::shared_ptr<Realizations>& R, std::shared_ptr<AbstractModel>& M,
                                     const std::shared_ptr<Data>& D, int IterationNumber) = 0;


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
    
    /// Candidates random variables, corresponding to those in the Model
    CandidateRandomVariables m_CandidateRandomVariables;
    


};


#endif //_AbstractSampler_h
