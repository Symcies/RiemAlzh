#ifndef _Algorithm_h
#define _Algorithm_h

#include <LongitudinalModel.h>
#include "../Samplers/AbstractSampler.h"
#include <chrono>

typedef std::vector<std::vector<std::pair<double, std::vector<double>>>> Data;


class Algorithm {
public:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    Algorithm();
    ~Algorithm();


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Getter(s) and Setter(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Get algorithm Parameters
    inline std::vector<double> GetParameters() const { return m_Model->GetAlgorithmParameters(); }

    // Set the longitudinal model used
    inline void SetModel(std::shared_ptr<LongitudinalModel> M) { m_Model = M; };

    // Set the data
    inline void SetData(Data *D) { m_Data = D; };

    // Set the sampler
    inline void SetSampler(std::shared_ptr<AbstractSampler> S) { m_Sampler = S; };


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Initialize all the attributes, in the Algo class as well as Longitudinal Model and Manifold
    void Initialize();

    // Simulation step of the MCMC SAEM
    void ComputeMCMCSAEM(int NumberOfIterations);

protected:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Longitudinal model
    std::shared_ptr<LongitudinalModel> m_Model;

    // Data
    Data *m_Data;

    // Sampler
    std::shared_ptr<AbstractSampler> m_Sampler;


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // MCMC SAEM Attribute(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Sufficient statistics of the MCMC SAEM algorithm S( y, z(k) )
    std::vector<std::vector<double>> m_SufficientStatistics;

    // Stochastic sufficient statistic approximation S_j(k)
    std::vector<std::vector<double>> m_StochasticSufficientStatistics;


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // MCMC SAEM Methods :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Compute the updates epsilon(k) to learn the stochastic sufficient statistics
    double DecreasingStepSize(int k, int NoMemory);

    // Compute the first step : simulation step
    void ComputeSimulationStep();

    // Compute the second step : Compute the sufficient statistics
    void ComputeSufficientStatistics();

    // Compute the third step : Compute the stochastic approximation step
    void ComputeStochasticApproximation(int );

    // Compute the fourth step : Maximization step
    void ComputeMaximizationStep();

};


#endif //_Algorithm_h
