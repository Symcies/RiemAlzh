#include "Algorithm.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

Algorithm
::Algorithm()
{
    m_Data = nullptr;
}

Algorithm
::~Algorithm()
{
    delete m_Data;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
Algorithm
::Initialize()
{
    // Initialize the Longitudinal model
    m_Model->SetData(m_Data);
    m_Model->Initialize();
    m_Model->Update();

    // Initialize stochastic sufficient statistics
    m_StochasticSufficientStatistics = m_Model->InitializeSufficientStochasticStatistics();
}

void
Algorithm
::ComputeMCMCSAEM(int NumberOfIterations) {

    for(int k = 0; k < NumberOfIterations ; ++k)
    {
        ComputeSimulationStep();
        ComputeSufficientStatistics();
        ComputeStochasticApproximation(k);
        ComputeMaximizationStep();
        GetParameters();
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// MCMC SAEM Methods :
////////////////////////////////////////////////////////////////////////////////////////////////////

double
Algorithm
::DecreasingStepSize(int k, int NoMemory)
{

    if(k<=NoMemory)
    {
        return 0;
    }
    else
    {
        return 1 / pow((k-NoMemory), 0.6);
    }
}

void
Algorithm
::ComputeSimulationStep()
{
    std::vector<RandomVariableToSample> RVToSample = m_Model->GetRandomVariableToSample();

    for(std::vector<RandomVariableToSample>::iterator it = RVToSample.begin() ; it != RVToSample.end() ; ++it)
    {
        bool q = m_Sampler->Sample(*it->first, *it->second, *m_Model);
        if(q)
        {
            m_Model->Update();
        }
    }

}

void
Algorithm
::ComputeSufficientStatistics()
{
    m_SufficientStatistics = m_Model->ComputeSufficientStatistics();
}

void
Algorithm
::ComputeStochasticApproximation(int k)
{
    double StepSize = DecreasingStepSize(k, 100);

    int i = 0;
    std::vector<std::vector< double >> NewStochasticStatistics;

    for(std::vector<std::vector<double>>::iterator it = m_StochasticSufficientStatistics.begin() ; it != m_StochasticSufficientStatistics.end() ; ++it)
    {
        std::vector<double> Coordinate;
        int j = 0;
        for(std::vector<double>::iterator it2 = it->begin() ; it2 != it->end() ; ++it2)
        {
            double Value = *it2 + StepSize * (m_SufficientStatistics[i][j] - *it2);
            Coordinate.push_back(Value);
            j += 1;
        }
        NewStochasticStatistics.push_back(Coordinate);
        i += 1;
    }

    m_StochasticSufficientStatistics = NewStochasticStatistics;

}

void
Algorithm
::ComputeMaximizationStep()
{
    m_Model->ComputeMaximizationStep(m_StochasticSufficientStatistics);
}



















