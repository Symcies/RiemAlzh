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

    /// Key : name of the parameter. Value : Value of the parameter
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
    GaussianRandomVariable GetRandomVariable(std::string rand_var_name, int num_real) const;
    GaussianRandomVariable GetRandomVariable(int rand_var_key, int num_real) const;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the candidate random variables
    void InitializeCandidateRandomVariables(const Realizations& reals, const AbstractModel& model);

    /// Update the variance of the random variable
    void UpdatePropositionVariableVariance(std::string name, int num_real, ScalarType new_variance);


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Proposition laws of the realizations
    PropositionDistribution proposition_distribution_;
    ProposDistrib new_proposition_distribution_;

    /// Key -> name conversion for the random variable names
    IntStringHash int_to_string_key_;

    /// name -> Key conversion for the random variable names
    StringIntHash string_to_int_key_;

};


#endif //_CandidateRandomVariables_h
