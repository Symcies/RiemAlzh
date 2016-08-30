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

    int NbMaxIterations = 50;
    InitializeRealization((int)D->size());
    InitializeStochasticSufficientStatistics(m_Model->GetSufficientStatistics(m_Realizations, D));

    for(int k = 0; k<NbMaxIterations; ++k)
    {
        std::cout << "--------------------- Iteration " << k << " -------------------------------" << std::endl;
        ComputeSimulationStep(D);
        std::vector< std::vector< double >> SufficientStatistics = m_Model->GetSufficientStatistics(m_Realizations, D);
        ComputeStochasticApproximation(k, SufficientStatistics);
        m_Model->UpdateRandomVariables(m_StochasticSufficientStatistics, D);
    }
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
    InitializeModelParameters(m_Realizations);
}

void
Algorithm
::InitializeModelParameters(std::shared_ptr<Realizations>& R)
{
    m_Model->InitializeModelParameters(R);
}

void
Algorithm
::ComputeSimulationStep(const std::shared_ptr<Data>& D)
{
    typedef Realizations::iterator ReaIter;
    typedef RandomVariableMap::iterator RandVarIter;

    m_Sampler->Sample(m_Realizations, m_Model, D);
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