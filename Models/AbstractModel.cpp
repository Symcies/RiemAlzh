#include "AbstractModel.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

AbstractManifold
::AbstractManifold()
{ }

AbstractManifold
::~AbstractManifold()
{ }


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


std::shared_ptr< AbstractRandomVariable > // RandomVariable
AbstractModel
::GetRandomVariable(std::string name)
{
    if(m_IndividualRandomVariables.count(name))
    {
        return m_IndividualRandomVariables.at(name);
    }
    else if(m_PopulationRandomVariables.count(name))
    {
        return m_PopulationRandomVariables.at(name);
    }
    else if(m_ManifoldRandomVariables.count(name))
    {
        return  m_ManifoldRandomVariables.at(name);
    }
    else
    {
        // TODO : change to a real warning
        std::cout << "OWN WARNING : The key does not exist in the map" << std::endl;
    }
}


std::map< std::string, std::shared_ptr< AbstractRandomVariable >>
AbstractModel
::GetRandomVariables()
{
    RandomVariableMap RandomVariables;

    RandomVariables.insert(m_ManifoldRandomVariables.begin(), m_ManifoldRandomVariables.end());
    RandomVariables.insert(m_PopulationRandomVariables.begin(), m_PopulationRandomVariables.end());
    RandomVariables.insert(m_IndividualRandomVariables.begin(), m_IndividualRandomVariables.end());

    return RandomVariables;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

std::map<std::string, std::vector<double>>
AbstractModel
::SimulateRealizations(int NumberOfSubjects)
{
    Realizations R;

    // Initialize the realization of the population-wide random variables.
    // They are shared among the individual thus sampled only once
    for(auto it : m_PopulationRandomVariables)
    {
        std::vector<double> Realization;
        Realization.push_back(it.second->Sample());
        R.insert(std::pair< std::string, std::vector<double>> (it.first, Realization));
    }

    // Initialize the realization of the individual random variables.
    // One realization is sampled for each individual, tagged with the i coefficient
    for(auto it : m_IndividualRandomVariables)
    {
        std::vector<double> Realization;
        for(int i = 0; i < NumberOfSubjects; ++i)
        {
            Realization.push_back(it.second->Sample());
        }
        R.insert(std::pair< std::string, std::vector<double>> (it.first, Realization));
    }

    // Initialize the realization of the manifold random variables/
    // They are shared among the individual thus sampled only once
    for(auto it : m_ManifoldRandomVariables)
    {
        std::vector<double> Realization;
        Realization.push_back(it.second->Sample());
        R.insert(std::pair< std::string, std::vector<double>> (it.first, Realization));
    }
    return R;
}