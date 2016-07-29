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


std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> // RandomVariable
AbstractModel
::GetRandomVariable(std::string name)
{
    if(m_IndividualRandomVariables.count(name))
    {
        RandomVariable R(name, m_IndividualRandomVariables.at(name));
        return R;
    }
    else if(m_PopulationRandomVariables.count(name))
    {
        RandomVariable R(name, m_PopulationRandomVariables.at(name) );
        return R;
    }
    else
    {
        std::cout << "The key does not exist in the map" << std::endl;
    }
}