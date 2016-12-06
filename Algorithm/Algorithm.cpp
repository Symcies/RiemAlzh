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

    int NbMaxIterations = 15000;
    InitializeModel(D);
    InitializeSampler();
    InitializeStochasticSufficientStatistics(m_Model->GetSufficientStatistics(m_Realizations, D));

//    double a1 = 0;
//    double a2 = 0;
//    double a3 = 0;
//    double a4 = 0;

    for(int k = 0; k<NbMaxIterations; ++k)
    {
        if( k%1 == 0 ) { std::cout  << std::endl << "--------------------- Iteration " << k << " -------------------------------" << std::endl; }
        //clock_t a = clock();
        //ComputeOutputs();
        ComputeSimulationStep(D, k);
        //clock_t b = clock();
        //ComputeOutputs();
        SufficientStatisticsVector SufficientStatistics = m_Model->GetSufficientStatistics(m_Realizations, D);
        //clock_t c = clock();
        //ComputeOutputs();
        ComputeStochasticApproximation(k, SufficientStatistics);
        //clock_t d = clock();
        //ComputeOutputs();
        m_Model->UpdateRandomVariables(m_StochasticSufficientStatistics, D);
        //clock_t e = clock();
        if( k%1 == 0 ) 
        { 
            ComputeOutputs();
            std::cout << "LogLikelihood : " << m_Model->ComputeLogLikelihood(m_Realizations, D) << std::endl; 
        }
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
    for(auto&& it : m_StochasticSufficientStatistics)
    {
        std::fill(it.begin(), it.end(), 0.0);
    }
}

void
Algorithm
::InitializeRealization(unsigned int NbIndividuals)
{

}


void
Algorithm
::InitializeModel(const std::shared_ptr<Data> D) 
{
    m_Model->Initialize(D);
    MultiRealizations R = m_Model->SimulateRealizations((int)D->size());
    m_Realizations = std::make_shared<MultiRealizations>(R);
    m_Model->UpdateParameters(m_Realizations);
    
    for(auto&& it : *m_Realizations)
    {
        m_OutputRealizations << it.first << ", ";

        VectorType v(it.second.size(), 0);
        m_AcceptanceRatios[it.first] = v;
    }
    m_OutputRealizations << std::endl;
}

void 
Algorithm
::InitializeSampler()
{
    m_Sampler->InitializeSampler(m_Realizations);
}

void
Algorithm
::ComputeSimulationStep(const std::shared_ptr<Data>& D, int Iteration)
{
    MultiRealizations&& R = m_Sampler->Sample(m_Realizations, m_Model, D, Iteration);
    ComputeAcceptanceRatio(R, Iteration);
    m_Realizations = std::make_shared<MultiRealizations>(R);
    
}


void 
Algorithm
::ComputeStochasticApproximation(double iteration, SufficientStatisticsVector SufficientStatistics)
{
    typedef std::vector< VectorType > SufficientStatisticsVector;
    SufficientStatisticsVector NewStochasticSufficientStatistics;

    double NoMemoryTime = 7500;  // TODO : Initialize, maybe out of the Compute function? Maybe in the decreasing step size function 
    double StepSize = DecreasingStepSize(iteration, NoMemoryTime);

    auto IterStat = SufficientStatistics.begin();
    auto IterStochStat = m_StochasticSufficientStatistics.begin();

    for(    ; IterStat != SufficientStatistics.end() && IterStochStat != m_StochasticSufficientStatistics.end()
            ; ++IterStat, ++IterStochStat)
    {
        VectorType S(IterStat->size());

        auto IterCoordStat = IterStat->begin();
        auto IterCoordStochStat = IterStochStat->begin();
        auto IterS = S.begin();
        for(    ; IterCoordStat != IterStat->end() && IterCoordStochStat != IterStochStat->end() && IterS != S.end()
                ; ++IterCoordStat, ++IterCoordStochStat, ++IterS)
        {
            *IterS = *IterCoordStochStat + StepSize * (*IterCoordStat - *IterCoordStochStat);
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
    /*
    for(const auto& it : *m_Realizations)
    {
        for(const auto& it2 : it.second)
        {
            m_OutputRealizations << it2 << ",";
        }
    }
    m_OutputRealizations << std::endl;
    */
    
    m_Model->ComputeOutputs();
}


void
Algorithm
::ComputeAcceptanceRatio(MultiRealizations& R, int Iteration)
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
          *IterAcceptRatio = (*IterAcceptRatio * Iteration + Change ) / (Iteration + 1);
      }
      
  }
    
    if(Iteration%1 == 0)
    {
        std::cout << "AcceptRatio: ";
        for(const auto& it : *m_Realizations)
        {
            std::cout << it.first << ": ";
            double Ave = 0;
            double Max = 0;
            double Min = 1;
            for(auto it2 : m_AcceptanceRatios.at(it.first))
            {
                Ave += it2;
                if(m_AcceptanceRatios.at(it.first).size() > 1)
                {
                    Max = std::max(it2, Max);
                    Min = std::min(it2, Min);
                }
            }
            std::cout << Ave/it.second.size() ;
            if(Min != 1) std::cout << " & " << Min;
            if(Max != 0) std::cout << " & " << Max;
            std::cout << ". ";
                                   
        }
        std::cout << std::endl;
    }
}



