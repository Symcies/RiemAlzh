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
const
{
    if(m_IndividualRandomVariables.count(name))
    {
        return m_IndividualRandomVariables.at(name);
    }
    else if(m_PopulationRandomVariables.count(name))
    {
        return m_PopulationRandomVariables.at(name);
    }
    else
    {
        // TODO : change to a real warning
        std::cout << "OWN WARNING : The key does not exist in the map" << std::endl;
    }
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////



Realizations
AbstractModel
::SimulateRealizations(int NumberOfSubjects)
{
    Realizations R;
    int Key = 0;

    // Initialize the realization of the population-wide random variables.
    // They are shared among the individual thus sampled only once
    for(auto it : m_PopulationRandomVariables)
    {
        VectorType Realization(1, it.second->Sample());
        R.AddRealizations(it.first, Key, Realization);
        ++Key;
    }

    // Initialize the realization of the individual random variables.
    // One realization is sampled for each individual, tagged with the i coefficient
    for(auto it : m_IndividualRandomVariables)
    {
        VectorType Realization(NumberOfSubjects);
        for(auto it2 = Realization.begin(); it2 != Realization.end(); ++it2)
        {
            double Sample = it.second->Sample(); 
            *it2 = Sample;
        }
        R.AddRealizations(it.first, Key, Realization);
        ++Key;
    }

    return R;
}