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
       m_Model->UpdateSubjectSpecificParameters(*m_Data);
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
        m_Model->UpdateSubjectSpecificParameters(*m_Data);
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
    // TO DO
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// MCMC SAEM Methods :
////////////////////////////////////////////////////////////////////////////////////////////////////

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
    m_SufficientStatistics.clear();

    // Compute S1 and S2
    std::vector<double> S1;
    std::vector<double> S2;

    int i = 0;
    for(Data::iterator it = m_Data->begin() ; it != m_Data->end() ; ++it)
    {
        for(std::vector<std::pair<double, std::vector<double>>>::iterator it2 = it->begin(); it2 != it->end(); ++it2)
        {
            std::vector<double> ParallelTransport = m_Model->ComputeParallelTransport(i, it2->first);
            double SumS1 = 0;
            int j = 0;
            for(std::vector<double>::iterator it3 = ParallelTransport.begin() ; it3 != ParallelTransport.end() ; ++it3)
            {
                SumS1 += *it3 * it2->second[j];
                j += 1;
            }
            S1.push_back(SumS1);
        }
        i += 1;
    }
    m_SufficientStatistics.push_back(S1);

    // Compute S2
    std::vector<double> S2;

}
























