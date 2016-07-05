#include "Algorithm.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

Algorithm
::Algorithm()
{
    m_InitiateModel = false;
}

Algorithm
::~Algorithm()
{
    delete m_Model;
    delete m_Data;
    delete m_Sampler;

}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Getter(s) and Setter(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
Algorithm
::SetModel(LongitudinalModel *M)
{
    m_Model = M;

    if(m_InitiateModel)
    {
       m_Model->UpdateSubjectSpecificParameters(m_Data);
    }
    else
    {
        m_InitiateModel = true;
    }
}

void
Algorithm
::SetData(Data *D)
{
    m_Data = D;

    if(m_InitiateModel)
    {
        m_Model->UpdateSubjectSpecificParameters(m_Data);
    }
    else
    {
        m_InitiateModel = true;
    }
}

void
Algorithm
::SetSampler(AbstractSampler *S)
{
    m_Sampler = S;
}

void
Algorithm
::ComputeMCMCSAEM() {
    // The parameters and random variables are updated automatically in the model

    // Need to initialize S
    m_Model->InitializeSufficientStochasticStatistics();

    // Epsilon(k) is calculated at each time step thanks to the DecreasingStepSize function
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
        m_Sampler->Sample(*it->first, *it->second, *m_Model);
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
::MaximizationStep()
{
    m_Model->ComputeMaximizationStep(m_StochasticSufficientStatistics);
}



















