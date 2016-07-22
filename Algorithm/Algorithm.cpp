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
        //ComputeSimulationStep();
        //Comp
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

    for(RandomVariableMap::iterator it = PopulationRandomVariable.begin() ; it != PopulationRandomVariable.end(); ++it)
    {
        m_Realization.insert( std::pair<std::string, double> (it->first, it->second->Sample()));
    }

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