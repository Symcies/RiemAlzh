#ifndef _MultiRandomVariables_h
#define _MultiRandomVariables_h

#include <unordered_map>
#include <memory>


#include "Realizations.h"
#include "AbstractRandomVariable.h"
#include "GaussianRandomVariable.h"

typedef double ScalarType;


class MultiRandomVariables {

public:
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
    typedef typename std::unordered_map<int, std::shared_ptr<AbstractRandomVariable>> IntRandomVariableHash;
    typedef typename std::unordered_map<std::string, int> StringIntHash;
    typedef typename std::unordered_map<std::string, ScalarType > StringScalarHash;
    typedef typename std::unordered_map<int, ScalarType > IntScalarHash;
    typedef typename std::unordered_map<int, int > IntIntHash;
    typedef typename std::unordered_map<int, std::string> IntStringHash;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
   
    MultiRandomVariables();
    ~MultiRandomVariables();
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    std::unique_ptr<AbstractRandomVariable> GetRandomVariable(std::string Name) const;
    
    std::unique_ptr<AbstractRandomVariable> GetRandomVariable(int Key) const;
    
    void Clear();
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
    void AddRandomVariable(std::string Name, std::string Type, const std::vector<double>& Parameters);
    
    /// Update a random variable based on its name
    void UpdateRandomVariable(std::string Name, IntScalarHash Parameters);
    void UpdateRandomVariable(std::string Name, StringScalarHash Parameters);
    
    /// Update a random variable based on its key
    void UpdateRandomVariable(int Key, IntScalarHash Parameters);
    void UpdateRandomVariable(int Key, StringScalarHash Parameters);
    
    /// Simulate realizations 
    Realizations SimulateRealizations(StringIntHash NumberOfRealizationsPerRandomVariable);
    
private:
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Hash table of the random variables
    IntRandomVariableHash m_RandomVariables;
    
    /// Convert the names into the int key
    StringIntHash m_StringToIntKey;
    
    /// Convert the int into their names
    IntStringHash m_IntToStringKey;
    
    /// Convert the variable key into their type
    IntStringHash m_KeyToRandomVariableStringType;
    IntIntHash    m_KeyToRandomVariableIntType;
    
    /// Int key counter
    int m_KeyCounter = 0;
    
};


#endif //_MultiRandomVariables_h
