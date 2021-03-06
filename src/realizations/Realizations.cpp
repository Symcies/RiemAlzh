#include "Realizations.h"


    
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////
    
Realizations
::Realizations() 
{
    
}

Realizations
::~Realizations() 
{
    
}




////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void
Realizations
::AddRealizations(std::string Name, int Key, VectorType& Values) 
{
    m_StringToIntKey.insert({Name, Key});
    m_IntToStringKey.insert({Key, Name});
    m_Realizations.insert({Key, Values});
}