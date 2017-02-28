#include "Algorithm.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

Algorithm
::Algorithm()
{
}

Algorithm
::Algorithm(io::AlgorithmSettings& Settings) 
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
::ComputeMCMCSAEM(const OldData& D, const Observations& Obs)
{
    InitializeModel(Obs);
    InitializeSampler();
    InitializeStochasticSufficientStatistics(Obs);

    for(m_IterationCounter = 0; m_IterationCounter < m_MaxNumberOfIterations; m_IterationCounter += 1)
    {
        if( m_IterationCounter%m_CounterToDisplayOutputs == 0 ) { std::cout  << std::endl << "--------------------- Iteration " << m_IterationCounter << " -------------------------------" << std::endl; }
        
        ComputeSimulationStep(D, Obs);
        SufficientStatisticsVector SufficientStatistics = m_Model->GetSufficientStatistics(*m_Realizations, Obs);
        ComputeStochasticApproximation(SufficientStatistics);
        m_Model->UpdateRandomVariables(m_StochasticSufficientStatistics);
        
        if( m_IterationCounter%m_CounterToDisplayOutputs == 0 ) { DisplayOutputs(); }
        if( m_IterationCounter%m_CounterToSaveData == 0) { m_Model->SaveData(m_IterationCounter, *m_Realizations); }
        
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
Algorithm
::InitializeStochasticSufficientStatistics(const Observations& Obs)
{
    m_StochasticSufficientStatistics = m_Model->GetSufficientStatistics(*m_Realizations, Obs);
    for(auto&& it : m_StochasticSufficientStatistics)
    {
        std::fill(it.begin(), it.end(), 0.0);
    }
}


void
Algorithm
::InitializeModel(const Observations& Obs) 
{
    
    m_Model->Initialize(Obs);
    Realizations R = m_Model->SimulateRealizations();
    
    m_Realizations = std::make_shared<Realizations>(R);
    m_Model->UpdateModel(R, -1);
    
    for(auto it = m_Realizations->begin(); it != m_Realizations->end(); ++it)
    {
        VectorType v(it->second.size(), 0);
        m_AcceptanceRatios[it->first] = v;
    }

}

void 
Algorithm
::InitializeSampler()
{
    m_Sampler->InitializeSampler(*m_Realizations, *m_Model);
}

void
Algorithm
::ComputeSimulationStep(const OldData& D, const Observations& Obs)
{
    Realizations PreviousRealisations2 = *m_Realizations;
    m_Sampler->Sample(*m_Realizations, *m_Model, D, Obs);
    ComputeAcceptanceRatio(PreviousRealisations2);
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
::ComputeAcceptanceRatio(Realizations& PreviousReals)
{ 
    

  for(auto it = m_Realizations->begin(); it != m_AcceptanceRatios.end(); ++it)
  {
      int KeyVariable = it->first;
      
      VectorType NewReal = it->second;
      VectorType PrevReal = PreviousReals.at(KeyVariable);
      
      auto IterPrevReal = PrevReal.begin();
      auto IterNewReal = NewReal.begin();
      auto IterAcceptRatio = m_AcceptanceRatios.at(KeyVariable).begin();
      
      for(    ; IterPrevReal != PrevReal.end() && IterNewReal != NewReal.end() && IterAcceptRatio != m_AcceptanceRatios.at(KeyVariable).end()
              ; ++IterPrevReal, ++IterNewReal, ++IterAcceptRatio)
      {
          bool Change = (*IterNewReal != *IterPrevReal);
          *IterAcceptRatio = (*IterAcceptRatio * m_IterationCounter + Change ) / (m_IterationCounter + 1);
      }
  }
    
    if(m_IterationCounter%m_CounterToDisplayOutputs == 0) { DisplayAcceptanceRatio(); }
}

void
Algorithm
::DisplayAcceptanceRatio() {
    std::cout << "AcceptRatio: ";
    
    auto NamesToShow = {"P", "Tau", "Ksi", "Beta#1", "Delta#3", "Nu"};
    
    for(auto it = NamesToShow.begin(); it != NamesToShow.end(); ++it)
    {
        std::string Name = *it;
        int Key = m_Realizations->ReverseNameToKey(Name);
        VectorType Ratios = m_AcceptanceRatios.at(Key);
        
        std::cout << Name << ": " << Ratios.mean_value();
        if(Ratios.size() != 1)
            std::cout << " & Min: " << Ratios.min_value() << " & Max: " << Ratios.max_value();
        std::cout << ". ";
    }
    std::cout << std::endl;
    
    /// Useless for now because all delta or beta are the same
    /*
    auto InName = {"Delta", "Beta"};
    for(auto it = InName.begin(); it != InName.end(); ++it)
    {
        double Min = 1;
        double Max = 0;
        double Mean = 0;
        int Count = 0;
        for(auto it2 = m_AcceptanceRatios.begin(); it2 != m_AcceptanceRatios.end(); ++it2)
        {
            std::string Name = m_Realizations->ReverseKeyToName(it2->first);
            Name = Name.substr(0, Name.find_first_of("#"));
            if(Name == *it)
            {
                ++Count;
                double AccepVal = it2->second(0);
                Min = std::min(Min, AccepVal);
                Max = std::max(Max, AccepVal);
                Mean += AccepVal;
            }
        }
        
        std::cout << *it <<"(Min/Mean/Max): " << Min << "/" << Mean/Count << "/" << Max << ". ";  
    }
    std::cout << std::endl;
    */
    
}

