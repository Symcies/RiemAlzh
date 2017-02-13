#include "MultiRandomVariables.h"




////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////
   
MultiRandomVariables
::MultiRandomVariables()
{
    
}

MultiRandomVariables
::~MultiRandomVariables() 
{
    
}
    
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<AbstractRandomVariable>
MultiRandomVariables
::GetRandomVariable(std::string Name) 
const 
{
    if(m_KeyToRandomVariableStringType.at(m_StringToIntKey.at(Name)) == "Gaussian") 
    {
        // TODO : Change to GetParameters(0) and GetPArameters(1) 
        ScalarType Mean = m_RandomVariables.at(m_StringToIntKey.at(Name))->GetParameter("Mean");
        ScalarType Variance = m_RandomVariables.at(m_StringToIntKey.at(Name))->GetParameter("Variance");
        return std::make_unique<GaussianRandomVariable>(Mean, Variance);
    }
}

std::unique_ptr<AbstractRandomVariable>
MultiRandomVariables
::GetRandomVariable(int Key) 
const 
{
    if(m_KeyToRandomVariableIntType.at(Key)) 
    {
        // TODO : Change to GetParameters(0) and GetPArameters(1) 
        ScalarType Mean = m_RandomVariables.at(Key)->GetParameter("Mean");
        ScalarType Variance = m_RandomVariables.at(Key)->GetParameter("Variance");
        return std::make_unique<GaussianRandomVariable>(Mean, Variance);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void
MultiRandomVariables
::AddRandomVariable(std::string Name, std::string Type, const std::vector<double>& Parameters) 
{
    std::shared_ptr<AbstractRandomVariable> RV;
    if(Type == "Gaussian") 
    {
        RV = std::make_shared<GaussianRandomVariable>(Parameters[0], Parameters[1]);
        m_KeyToRandomVariableStringType.insert({m_KeyCounter, "Gaussian"});
        m_KeyToRandomVariableIntType.insert({m_KeyCounter, 0});
    }
    
    m_StringToIntKey.insert({Name, m_KeyCounter});
    m_IntToStringKey.insert({m_KeyCounter, Name});
    m_RandomVariables.insert({m_KeyCounter, RV});
    ++m_KeyCounter;
}

void 
MultiRandomVariables
::UpdateRandomVariable(std::string Name, IntScalarHash Parameters) 
{
    if(m_RandomVariables.find(m_StringToIntKey.at(Name)) != m_RandomVariables.end())
        m_RandomVariables.at(m_StringToIntKey.at(Name))->Update(Parameters);
    else
        std::cerr << Name << " is not found in the random variables";
}

void
MultiRandomVariables
::UpdateRandomVariable(std::string Name, StringScalarHash Parameters) 
{
    if(m_RandomVariables.find(m_StringToIntKey.at(Name)) != m_RandomVariables.end())
        m_RandomVariables.at(m_StringToIntKey.at(Name))->Update(Parameters);
    else
        std::cerr << Name << " is not found in the random variables";
}


void
MultiRandomVariables
::UpdateRandomVariable(int Key, IntScalarHash Parameters) 
{
    if(m_RandomVariables.find(Key) != m_RandomVariables.end())
        m_RandomVariables.at(Key)->Update(Parameters);
    else
        std::cerr << m_IntToStringKey.at(Key) << " is not found in the random variable";
}


void 
MultiRandomVariables
::UpdateRandomVariable(int Key, StringScalarHash Parameters) 
{
    if(m_RandomVariables.find(Key) != m_RandomVariables.end())
        m_RandomVariables.at(Key)->Update(Parameters);
    else
        std::cerr << m_IntToStringKey.at(Key) << " is not found in the random variable";
}


Realizations
MultiRandomVariables
::SimulateRealizations(StringIntHash NumberOfRealizationsPerRandomVariable) 
{
    Realizations R;
    
    for(auto it = NumberOfRealizationsPerRandomVariable.begin(); it != NumberOfRealizationsPerRandomVariable.end(); ++it)
    {
        std::string Name = it->first;
        int Key = m_StringToIntKey.at(Name);
        int NumberOfRealizations = it->second;
        
        VectorType Real = m_RandomVariables.at(Key)->Samples(NumberOfRealizations);
        R.AddRealizations(Name, Key, Real);
    }
    
    return R;
}