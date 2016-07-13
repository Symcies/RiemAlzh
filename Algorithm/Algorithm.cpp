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

    double SimulationTime = 0;
    double SuffStatTime = 0;
    double StochApproxTime = 0;
    double MaxStepTime = 0;

    for(int k = 0; k < NumberOfIterations ; ++k)
    {
        auto t1 = std::chrono::system_clock::now();
        ComputeSimulationStep();
        auto t2 = std::chrono::system_clock::now();
        ComputeSufficientStatistics();
        auto t3 = std::chrono::system_clock::now();
        ComputeStochasticApproximation(k);
        auto t4 = std::chrono::system_clock::now();
        ComputeMaximizationStep();
        auto t5 = std::chrono::system_clock::now();
        //GetParameters();

        auto p1 = t2 - t1;
        auto p2 = t3 - t2;
        auto p3 = t4 - t3;
        auto p4 = t5 - t4;
        SimulationTime += std::chrono::duration_cast<std::chrono::nanoseconds>(p1).count();
        SuffStatTime += std::chrono::duration_cast<std::chrono::nanoseconds>(p2).count();
        StochApproxTime += std::chrono::duration_cast<std::chrono::nanoseconds>(p3).count();
        MaxStepTime += std::chrono::duration_cast<std::chrono::nanoseconds>(p4).count();
    }

    std::cout << "Simulation Time : " << SimulationTime/1000 << std::endl;
    std::cout << "Sufficient Time : " << SuffStatTime/1000 << std::endl;
    std::cout << "Stochastic Time : " << StochApproxTime/1000 << std::endl;
    std::cout << "Maximizati Time : " << MaxStepTime/1000 << std::endl;
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
        std::cout << "-----------" << q << "-------------------" << std::endl;
    }
    std::cout << "----------------------------------------------------------------------------------------" << std::endl;

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
    double StepSize = DecreasingStepSize(k, 50);

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



















