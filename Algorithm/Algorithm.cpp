#include "Algorithm.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

Algorithm
::Algorithm()
{
    m_OutputRealizations.open("Realizations.txt",  std::ofstream::out | std::ofstream::trunc );
}

Algorithm
::~Algorithm()
{ }


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
Algorithm
::ComputeMCMCSAEM(const std::shared_ptr<Data>& D)
{

    int NbMaxIterations = 300;
    InitializeRealization((int)D->size());
    InitializeModel(m_Realizations);
    InitializeCandidateRandomVariables();
    InitializeStochasticSufficientStatistics(m_Model->GetSufficientStatistics(m_Realizations, D));

//    double a1 = 0;
//    double a2 = 0;
//    double a3 = 0;
//    double a4 = 0;

    for(int k = 0; k<NbMaxIterations; ++k)
    {
        std::cout  << std::endl << "--------------------- Iteration " << k << " -------------------------------" << std::endl;
        //clock_t a = clock();
        ComputeSimulationStep(D, k);
        //clock_t b = clock();
        std::vector< std::vector< double >> SufficientStatistics = m_Model->GetSufficientStatistics(m_Realizations, D);
        //clock_t c = clock();
        ComputeStochasticApproximation(k, SufficientStatistics);
        //clock_t d = clock();
        m_Model->UpdateRandomVariables(m_StochasticSufficientStatistics, D);
        //clock_t e = clock();
        ComputeOutputs();
//        a1 += b - a;
//        a2 += c - b;
//        a3 += d - c;
//        a4 += e - d;
    }

//    std::cout << "Simulate : " << double(a1) << std::endl;
//    std::cout << "SuffStat : " << double(a2) << std::endl;
//    std::cout << "Stochast : " << double(a3) << std::endl;
//    std::cout << "Maximize : " << double(a4) << std::endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
Algorithm
::InitializeStochasticSufficientStatistics(const SufficientStatisticsVector& S)
{
    m_StochasticSufficientStatistics = S;
    for(auto it = m_StochasticSufficientStatistics.begin(); it != m_StochasticSufficientStatistics.end(); ++it)
    {
        std::fill(it->begin(), it->end(), 0.0);
    }
}

void
Algorithm
::InitializeRealization(unsigned int NbIndividuals)
{
    Realizations R = m_Model->SimulateRealizations(NbIndividuals);
    m_Realizations = std::make_shared<Realizations>(R);

    for(auto it = m_Realizations->begin(); it != m_Realizations->end(); ++it)
    {
        m_OutputRealizations << it->first << ", ";

        std::vector<double> v(it->second.size(), 0);
        m_AcceptanceRatios[it->first] = v;
    }
    m_OutputRealizations << std::endl;
}

void
Algorithm
::InitializeCandidateRandomVariables()
{
    m_CandidateRandomVariables = std::make_shared<CandidateRandomVariables>();
    m_CandidateRandomVariables->InitializeCandidateRandomVariables(m_Model);
}

void
Algorithm
::InitializeModel(std::shared_ptr<Realizations> &R) 
{
    m_Model->UpdateParameters(R);
}

void
Algorithm
::ComputeSimulationStep(const std::shared_ptr<Data>& D, int Iteration)
{
    typedef Realizations::iterator ReaIter;
    typedef RandomVariableMap::iterator RandVarIter;

    Realizations R = m_Sampler->Sample(m_Realizations, m_Model, m_CandidateRandomVariables, D);
    ComputeAcceptanceRatio(R, Iteration);
    m_Realizations = std::make_shared<Realizations>(R);
    
}


void 
Algorithm
::ComputeStochasticApproximation(double iteration, std::vector< std::vector< double >> SufficientStatistics)
{
    typedef std::vector< std::vector< double >> SufficientStatisticsVector;
    SufficientStatisticsVector NewStochasticSufficientStatistics;

    double NoMemoryTime = 600;  // TODO : Initialize, maybe out of the Compute function? Maybe in the decreasing step size function 
    double StepSize = DecreasingStepSize(iteration, NoMemoryTime);

    auto IterStat = SufficientStatistics.begin();
    auto IterStochStat = m_StochasticSufficientStatistics.begin();

    for(    ; IterStat != SufficientStatistics.end() && IterStochStat != m_StochasticSufficientStatistics.end()
            ; ++IterStat, ++IterStochStat)
    {
        std::vector<double> S;

        auto IterCoordStat = IterStat->begin();
        auto IterCoordStochStat = IterStochStat->begin();
        for(    ; IterCoordStat != IterStat->end() && IterCoordStochStat != IterStochStat->end()
                ; ++IterCoordStat, ++IterCoordStochStat)
        {
            S.push_back( *IterCoordStochStat + StepSize * (*IterCoordStat - *IterCoordStochStat) );
        }

        NewStochasticSufficientStatistics.push_back(S);

    }
    
    m_StochasticSufficientStatistics = NewStochasticSufficientStatistics;

}


double
Algorithm
::DecreasingStepSize(double Iteration, double NoMemoryTime)
{
    double Epsilon = std::max(1.0, Iteration - NoMemoryTime);

    return 1.0 / pow(Epsilon, 0.6); // TODO : TO CHECK
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Output(s)
////////////////////////////////////////////////////////////////////////////////////////////////////


void
Algorithm
::ComputeOutputs()
{
    for(auto it : *m_Realizations)
    {
        for(auto it2 : it.second)
        {
            m_OutputRealizations << it2 << ",";
        }
    }
    m_OutputRealizations << std::endl;

    
}


void
Algorithm
::ComputeAcceptanceRatio(Realizations& R, int Iteration)
{
    std::cout << "AcceptRatio: ";
  for(auto it : *m_Realizations)
  {
      std::string NameVariable = it.first;
      std::cout << it.first << ": ";
      double AverageAcceptanceRatioToPrint = 0;
      
      std::vector<double> PrevReal = it.second;
      std::vector<double> NewReal = R.at(NameVariable);
      
      auto IterPrevReal = PrevReal.begin();
      auto IterNewReal = NewReal.begin();
      auto IterAcceptRatio = m_AcceptanceRatios.at(NameVariable).begin();
      
      for(    ; IterPrevReal != PrevReal.end() && IterNewReal != NewReal.end() && IterAcceptRatio != m_AcceptanceRatios.at(NameVariable).end()
              ; ++IterPrevReal, ++IterNewReal, ++IterAcceptRatio)
      {
          //std::cout << "b/a : " << *IterPrevReal << "/" << *IterNewReal << "->";
          bool Change = (*IterNewReal != *IterPrevReal);
          *IterAcceptRatio = (*IterAcceptRatio * Iteration + Change ) / (Iteration + 1);
          AverageAcceptanceRatioToPrint += *IterAcceptRatio;
      }
      
      std::cout << AverageAcceptanceRatioToPrint / PrevReal.size() << ". ";
  }
    std::cout << std::endl;
    
}



