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
::ComputeMCMCSAEM(std::shared_ptr<Data> D)
{
    int NbMaxIterations = 100;
    InitializeRealization(D->size());

    for(int k = 0; k<NbMaxIterations; ++k)
    {
        ComputeSimulationStep();
        std::vector< std::vector< double >> SufficientStatistics = m_Model->GetSufficientStatistics(m_Realization, *D);
        ComputeStochasticApproximation(k, SufficientStatistics);
        m_Model->UpdateRandomVariables(m_StochasticSufficientStatistics, *D);
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
Algorithm
::InitializeRealization(unsigned int NbIndividuals)
{
    m_Realization.clear();

    RandomVariableMap PopulationRandomVariable = m_Model->GetPopulationRandomVariables();
    RandomVariableMap IndividualRandomVariable = m_Model->GetIndividualRandomVariables();

    // Initialize the realization of the population-wide random variables.
    // They are shared among the individual thus sampled only once
    for(RandomVariableMap::iterator it = PopulationRandomVariable.begin() ; it != PopulationRandomVariable.end(); ++it)
    {
        m_Realization.insert( std::pair<std::string, double> (it->first, it->second->Sample()));
    }

    // Initialize the realization of the individual random variables.
    // One realization is sampled for each individual, tagged with the i coefficient
    for(RandomVariableMap::iterator it = IndividualRandomVariable.begin() ; it != IndividualRandomVariable.end() ; ++it)
    {
        for(int i = 0; i < NbIndividuals ; ++i)
        {
            m_Realization.insert( std::pair<std::string, double> (it->first + std::to_string(i), it->second->Sample()));
        }
    }

}

void
Algorithm
::ComputeSimulationStep()
{
    typedef Realizations::iterator ReaIter;
    typedef RandomVariableMap::iterator RandVarIter;

    for(Realizations::iterator it = m_Realization.begin(); it != m_Realization.end() ; ++it)
    {
        std::string key = it->first;
        double value = it->second;

        //m_Sampler->

    }

    // For each random variables, get the corresponding realization (double)
    // To each random variables, associate a proposition random variables (Maybe in the model?)
    // Get a realization according to the proposition variable

    // Loop over all the variables to sample
    // Get their realization
    // Get the associated proposition random variables
    // Get their realization as well
        // m_Sampler->Sampler(CurrentRandomVar, CurrentRealization, CandidateRandomVar, CandidateRealization);

}



void 
Algorithm
::ComputeStochasticApproximation(double iteration, std::vector< std::vector< double >> SufficientStatistics)
{
    typedef std::vector< std::vector< double >> SufficientStatisticsVector;
    SufficientStatisticsVector NewStochasticSufficientStatistics;

    double NoMemoryTime = 100;  // TODO : Initialize, maybe out of the Compute function? Maybe in the decreasing step size function 
    double StepSize = DecreasingStepSize(iteration, 100);

    for(std::pair<SufficientStatisticsVector::iterator, SufficientStatisticsVector::iterator> i(SufficientStatistics.begin(), m_StochasticSufficientStatistics.begin()) ;
        i.first != SufficientStatistics.end() && i.second != m_StochasticSufficientStatistics.end() ;
        ++i.first, ++i.second)
    {
        std::vector<double> S;

        for(std::pair<std::vector<double>::iterator, std::vector<double>::iterator> j(i.first->begin(), i.second->begin()) ;
            j.first != i.first->end() && j.second != i.second->end() ;
            ++j.first, ++j.second)
        {
            S.push_back(*j.second + StepSize * (*j.first - *j.second)); // TODO : CHECK WHAT i AND j are
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