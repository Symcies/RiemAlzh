#ifndef _Algorithm_h
#define _Algorithm_h

#include <iostream>
#include <fstream>
#include <cassert>


#include "AbstractModel.h"
#include "AbstractSampler.h"
#include "CandidateRandomVariables.h"

class Algorithm {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    typedef typename LinearAlgebra<ScalarType>::MatrixType MatrixType;
    typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
    
    typedef std::vector< std::vector< std::pair< VectorType, double> > > Data;
    typedef std::map<std::string, VectorType> Realizations;
    typedef std::vector<VectorType> SufficientStatisticsVector;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    Algorithm();
    
    Algorithm(unsigned int MaxNumberOfIterations, unsigned int BurnIn);
    
    ~Algorithm();



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    inline void SetModel(const std::shared_ptr<AbstractModel>& M ) { m_Model = M; }

    inline void SetSampler(std::shared_ptr<AbstractSampler> S) { m_Sampler = S; }

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compute the MCMC SAEM algorithm
    void ComputeMCMCSAEM(const Data& D);


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Method(s) of the MCMC SAEM:
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the stochastic approximation
    void InitializeStochasticSufficientStatistics(const Data& D);
        
    /// Initialize sampler
    void InitializeSampler();
    
    /// Initialize Manifold
    void InitializeModel(const Data& D);

    /// Compute the simulation step : Gibbs Sampling
    void ComputeSimulationStep(const Data& D);

    /// Compute the stochastic coefficient 
    void ComputeStochasticApproximation(SufficientStatisticsVector& S);

    /// Compute the decreasing step size of the approximation step
    double DecreasingStepSize();


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Output(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Number of iterations to wait till next output display
    unsigned int m_CounterToDisplayOutputs = 100;
    
    /// Number of iterations to wait till next data saving
    unsigned int m_CounterToSaveData = 200;
    
    /// Acceptance Ratios
    std::map<std::string, VectorType> m_AcceptanceRatios;
    
    /// Compute the acceptance ratio for each random variable
    void ComputeAcceptanceRatio(Realizations& R);
    
    /// Display acceptance ratio
    void DisplayAcceptanceRatio(Realizations& R);
    
    /// Display Outputs
    void DisplayOutputs();
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Abstract Model
    std::shared_ptr<AbstractModel> m_Model;

    /// Realisation of the random variables of the model
    std::shared_ptr<Realizations> m_Realizations;

    /// Abstract Sampler - for the MCMC SAEM 
    std::shared_ptr<AbstractSampler> m_Sampler;

    /// Stochastic Sufficient Statistics used in the stochastic approximation step
    SufficientStatisticsVector m_StochasticSufficientStatistics;
    
    /// Total number of iterations
    unsigned int m_MaxNumberOfIterations = 5001;
    
    /// Number of burn-in iterations
    unsigned int m_BurnIn = 10002;
    
    /// Number of iterations done by the MCMC-SAEM
    unsigned int m_IterationCounter = 0;
};


#endif //_Algorithm_h
