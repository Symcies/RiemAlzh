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
        std::vector< std::vector< double >> SufficientStatistics = m_Model->GetSufficientStatistics(&R, &D);
        ComputeStochasticSufficientCoefficient(k, SufficientStatistics);
        ComputeMaximizationStep();
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

    RandomVariableMap PopulationRandomVariable = m_Model->GetPopulationrandomVariables();
    RandomVariableMap IndividualRandomVariable = m_Model->GetIndividualRandomVariables();

    // Initialize the realization of the population-wide random variables.
    // They are shared among the individual thus sampled only once
    for(RandomVariableMap::iterator it = PopulationRandomVariable.begin() ; it != PopulationRandomVariable.end(); ++it)
    {
        m_Realization.insert( std::pair<std::string, double> (it->first, it->second->Sample()));
    }

    // Initialize the realization of the individual random variables.
    // One realization is sampled for each individual, tagged with the i coefficent
    for(RandomVariableMap::iterator it = IndividualRandomVariable.begin() ; it != IndividualRandomVariable.end() ; ++it)
    {
        for(int i = 0; i < NbIndividuals ; ++i)
        {
            std::string name = it->first. + i;
            name += i;
            m_Realization.insert( std::pair<std::string, double> (name, it->second->Sample()));
        }
    }

}

void
Algorithm
::ComputeSimulationStep()
{

    
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
    typedef std::vector< std::vector< double >> SuffStat;
    SuffStats NewStochasticSufficientStatistics;

    double NoMemoryTime = 100;  // TODO : Initialize, maybe out of the Compute function? Maybe in the decreasing step size function 
    StepSize = DecreasingStepSize(k, 100);

    for(std::pair<SuffStats, SuffStats> i(SuffficientStatistics.begin(), m_StochasticSufficientStatistics.begin()) ;
        i.first != SufficientStatistics.end() && i.second != m_StochasticSufficientStatistics.end() ;
        ++i.first, ++i.second)
    {
        std::vector<double> S;

        for(std::pair<std::vector<double>, std::vector<double> j(i.first.begin(), i.second.begin()) ;
            j.first != i.first.end() && j.second != i.first.end() ;
            ++j.first, ++j.second)
        {
            S.push_back(*j + StepSize * (*i - *j));
        }

        NewStochasticSufficientStatistics.push_back(S)
    }

    m_StochasticSufficientStatistics = NewStochasticSufficientStatistics
}



void
Algorithm
::DecreasingStepSize(double Iteration, double NoMemoryTime)
{
    double epsilon = max(1, Iteration - NoMemoryTime);

    return 1.0 / pow(Epsilon, 0.6) // TODO : TO CHECK 
}