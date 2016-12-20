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
    else
    {
        // TODO : change to a real warning
        std::cout << "OWN WARNING : The key does not exist in the map" << std::endl;
    }
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////



AbstractModel::Realizations
AbstractModel
::SimulateRealizations(int NumberOfSubjects)
{
    Realizations R;

    // Initialize the realization of the population-wide random variables.
    // They are shared among the individual thus sampled only once
    for(auto it : m_PopulationRandomVariables)
    {
        VectorType Realization(1, it.second->Sample());
        R.insert(std::pair< std::string, VectorType> (it.first, Realization));
        if(it.first.substr(0, it.first.find_first_of("#")) != "Beta")
        {
            double a = 0;
        }
    }

    // Initialize the realization of the individual random variables.
    // One realization is sampled for each individual, tagged with the i coefficient
    for(auto it : m_IndividualRandomVariables)
    {
        VectorType Realization(NumberOfSubjects);
        double Var = 0.0, Mean = 0.0;
        for(auto it2 = Realization.begin(); it2 != Realization.end(); ++it2)
        {
            double Sample = it.second->Sample(); 
            *it2 = Sample;
            Mean += Sample;
            Var += Sample*Sample ;
        }
        if(it.first == "Ksi")
        {
            std::cout << "Ksi Mean : " << Mean / NumberOfSubjects << std::endl;
            std::cout << "Ksi Var : " <<  (Var - Mean * Mean /NumberOfSubjects)/NumberOfSubjects << std::endl; 
        }
        if(it.first == "Tau")
        {
            std::cout << "Tau Mean : " << Mean / NumberOfSubjects << std::endl;
            std::cout << "Tau Var : " << (Var - Mean * Mean /NumberOfSubjects) / NumberOfSubjects << std::endl;
        }
        R.insert(std::pair< std::string, VectorType> (it.first, Realization));
    }

    return R;
}