#include "Algorithm.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

Algorithm
::Algorithm()
{ }

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

    int NbMaxIterations = 10;
    InitializeRealization((int)D->size());
    InitializeCandidateRandomVariables();
    InitializeStochasticSufficientStatistics(m_Model->GetSufficientStatistics(m_Realizations, D));

    double a1 = 0;
    double a2 = 0;
    double a3 = 0;
    double a4 = 0;

    for(int k = 0; k<NbMaxIterations; ++k)
    {
        std::cout  << std::endl << "--------------------- Iteration " << k << " -------------------------------" << std::endl;
        //clock_t a = clock();
        ComputeSimulationStep(D);
        //clock_t b = clock();
        std::vector< std::vector< double >> SufficientStatistics = m_Model->GetSufficientStatistics(m_Realizations, D);
        //clock_t c = clock();
        ComputeStochasticApproximation(k, SufficientStatistics);
        //clock_t d = clock();
        m_Model->UpdateRandomVariables(m_StochasticSufficientStatistics, D);
        //clock_t e = clock();
        ComputeRealizationsEvolution();
        //a1 += b - a;
        //a2 += c - b;
        //a3 += d - c;
        //a4 += e - d;
    }

    //std::cout << "Simulate : " << double(a1) << std::endl;
    //std::cout << "SuffStat : " << double(a2) << std::endl;
    //std::cout << "Stochast : " << double(a3) << std::endl;
    //std::cout << "Maximize : " << double(a4) << std::endl;
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
        m_RealizationsEvolution[it->first] = std::vector<std::vector<double>>(1, it->second);
    }
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
::ComputeSimulationStep(const std::shared_ptr<Data>& D)
{
    typedef Realizations::iterator ReaIter;
    typedef RandomVariableMap::iterator RandVarIter;

    m_Sampler->Sample(m_Realizations, m_Model, m_CandidateRandomVariables, D);
}


void 
Algorithm
::ComputeStochasticApproximation(double iteration, std::vector< std::vector< double >> SufficientStatistics)
{
    typedef std::vector< std::vector< double >> SufficientStatisticsVector;
    SufficientStatisticsVector NewStochasticSufficientStatistics;

    double NoMemoryTime = 100;  // TODO : Initialize, maybe out of the Compute function? Maybe in the decreasing step size function 
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
            if(isnan(*IterCoordStochStat + StepSize * (*IterCoordStat - *IterCoordStochStat)))
            {
                //std::cout << "Comput : " << *IterCoordStochStat << " & " << StepSize << " & " << *IterCoordStat << std::endl;
            }
            S.push_back( *IterCoordStochStat + StepSize * (*IterCoordStat - *IterCoordStochStat) );
        }

        NewStochasticSufficientStatistics.push_back(S);

    }

    //std::cout << "2. Stoch approx : " << NewStochasticSufficientStatistics[0][0] << " & " << NewStochasticSufficientStatistics[0][1] << std::endl;


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
::ComputeRealizationsEvolution()
{
    for(auto it = m_Realizations->begin(); it != m_Realizations->end(); ++it)
    {
        std::string Name = it->first;
        std::vector<double> Realization = it->second;
        m_RealizationsEvolution[Name].push_back(Realization);
    }
}