#include "Algorithm.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

Algorithm
::Algorithm()
{
}

Algorithm
::Algorithm(AlgorithmSettings& Settings) 
{
    m_MaxNumberOfIterations = Settings.GetMaximumNumberOfIterations();
    m_BurnIn = Settings.GetNumberOfBurnInIterations();
    m_CounterToDisplayOutputs = Settings.GetCounterToDisplayOutputs();
    m_CounterToSaveData = Settings.GetCounterToSaveData();
    
}


Algorithm
::Algorithm(unsigned int MaxNumberOfIterations, unsigned int BurnIn) 
{
    m_MaxNumberOfIterations = MaxNumberOfIterations;
    m_BurnIn = BurnIn;
}

Algorithm
::~Algorithm()
{ }


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
Algorithm
::ComputeMCMCSAEM(const Data& D)
{
    InitializeModel(D);
    InitializeSampler(D);
    InitializeStochasticSufficientStatistics(D);

    for(m_IterationCounter = 0; m_IterationCounter < m_MaxNumberOfIterations; m_IterationCounter += 1)
    {
        if( m_IterationCounter%m_CounterToDisplayOutputs == 0 ) { std::cout  << std::endl << "--------------------- Iteration " << m_IterationCounter << " -------------------------------" << std::endl; }
        
        ComputeSimulationStep(D);
        SufficientStatisticsVector SufficientStatistics = m_Model->GetSufficientStatistics(*m_Realizations, D);
        ComputeStochasticApproximation(SufficientStatistics);
        m_Model->UpdateRandomVariables(m_StochasticSufficientStatistics, D);
        
        if( m_IterationCounter%m_CounterToDisplayOutputs == 0 ) { DisplayOutputs(); }
        if( m_IterationCounter%m_CounterToSaveData == 0) { m_Model->SaveData(m_IterationCounter, *m_Realizations); }
        
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
Algorithm
::InitializeStochasticSufficientStatistics(const Data& D)
{
    m_StochasticSufficientStatistics = m_Model->GetSufficientStatistics(*m_Realizations, D);
    for(auto&& it : m_StochasticSufficientStatistics)
    {
        std::fill(it.begin(), it.end(), 0.0);
    }
}


void
Algorithm
::InitializeModel(const Data& D) 
{
    
    m_Model->Initialize(D);
    Realizations R = m_Model->SimulateRealizations((int)D.size());
    m_Realizations = std::make_shared<Realizations>(R);
    m_Model->UpdateModel(R, -1);
    
    for(auto&& it : *m_Realizations)
    {
        VectorType v(it.second.size(), 0);
        m_AcceptanceRatios[it.first] = v;
    }
}

void 
Algorithm
::InitializeSampler(const Data& D)
{
    m_Sampler->InitializeSampler(*m_Realizations, *m_Model, D);
}

void
Algorithm
::ComputeSimulationStep(const Data& D)
{
    Realizations RRR = *m_Realizations;
    Realizations&& R = m_Sampler->Sample(RRR, *m_Model, D, m_IterationCounter);
    ComputeAcceptanceRatio(R);
    m_Realizations = std::make_shared<Realizations>(R);
    
}


void 
Algorithm
::ComputeStochasticApproximation(SufficientStatisticsVector& S)
{   
    
    assert(S.size() == m_StochasticSufficientStatistics.size()); 
    
    double StepSize = DecreasingStepSize();
    auto itStochS = m_StochasticSufficientStatistics.begin();
    
    for(auto itS = S.begin(); itS != S.end(); ++itS, ++itStochS)
    {
        *itStochS += StepSize * (*itS - *itStochS);
    }
}


double
Algorithm
::DecreasingStepSize()
{
    double Q = (double)m_IterationCounter - (double)m_BurnIn;
    double Epsilon = std::max(1.0, Q);
    return 1.0 / pow(Epsilon, 0.6); // TODO : TO CHECK
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Output(s)
////////////////////////////////////////////////////////////////////////////////////////////////////


void
Algorithm
::DisplayOutputs()
{
    m_Model->DisplayOutputs(*m_Realizations);
}


void
Algorithm
::ComputeAcceptanceRatio(Realizations& R)
{ 
    

  for(const auto& it : *m_Realizations)
  {
      std::string NameVariable = it.first;
      
      VectorType PrevReal = it.second;
      VectorType NewReal = R.at(NameVariable);
      
      auto IterPrevReal = PrevReal.begin();
      auto IterNewReal = NewReal.begin();
      auto IterAcceptRatio = m_AcceptanceRatios.at(NameVariable).begin();
      
      for(    ; IterPrevReal != PrevReal.end() && IterNewReal != NewReal.end() && IterAcceptRatio != m_AcceptanceRatios.at(NameVariable).end()
              ; ++IterPrevReal, ++IterNewReal, ++IterAcceptRatio)
      {
          bool Change = (*IterNewReal != *IterPrevReal);
          *IterAcceptRatio = (*IterAcceptRatio * m_IterationCounter + Change ) / (m_IterationCounter + 1);
      }
      
  }
    
    if(m_IterationCounter%m_CounterToDisplayOutputs == 0) { DisplayAcceptanceRatio(R); }
}

void
Algorithm
::DisplayAcceptanceRatio(Realizations &R) {
    std::cout << "AcceptRatio: ";
    for (const auto &it : *m_Realizations) {
        std::cout << it.first << ": ";
        double Ave = 0;
        double Max = 0;
        double Min = 1;
        for (auto it2 : m_AcceptanceRatios.at(it.first)) {
            Ave += it2;
            if (m_AcceptanceRatios.at(it.first).size() > 1) {
                Max = std::max(it2, Max);
                Min = std::min(it2, Min);
            }
        }
        std::cout << Ave / it.second.size();
        if (Min != 1) std::cout << " & " << Min;
        if (Max != 0) std::cout << " & " << Max;
        std::cout << ". ";

    }
    std::cout << std::endl;
}

