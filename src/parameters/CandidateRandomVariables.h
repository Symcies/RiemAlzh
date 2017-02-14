#ifndef _CandidateRandomVariables_h
#define _CandidateRandomVariables_h

#include <map>

#include "LinearAlgebra.h"
#include "GaussianRandomVariable.h"
#include "AbstractRandomVariable.h"
#include "AbstractModel.h"
#include "Realizations.h"

class CandidateRandomVariables {

public:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
    
    /// Key : Name of the parameter. Value : Value of the parameter
    typedef std::map< std::string, std::vector< GaussianRandomVariable >> PropositionDistribution;
    typedef std::unordered_map< int, std::vector<GaussianRandomVariable>> ProposDistrib;
    typedef typename std::unordered_map<std::string, int> StringIntHash;
    typedef typename std::unordered_map<int, std::string> IntStringHash;


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    CandidateRandomVariables(); // TODO : Add one argument : the XML file
    ~CandidateRandomVariables();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Get a candidate random variable
    GaussianRandomVariable GetRandomVariable(std::string NameRandomVariable, int RealizationNumber) const;
    GaussianRandomVariable GetRandomVariable(int RandomVariableKey, int RealizationNumber) const;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the candidate random variables
    void InitializeCandidateRandomVariables(const Realizations& R);



protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Return the initial proposition distribution corresponding to the variable
    GaussianRandomVariable ReadInitialPropositionDistribution(std::string NameRandomVariable, ScalarType CurrentState);
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Proposition laws of the realizations
    PropositionDistribution m_PropositionDistribution;
    ProposDistrib m_NewPropositionDistribution;
    
    /// Key -> Name conversion for the random variable names
    IntStringHash m_IntToStringKey;
    
    /// Name -> Key conversion for the random variable names
    StringIntHash m_StringToIntKey;
    
};


#endif //_CandidateRandomVariables_h
