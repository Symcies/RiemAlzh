#ifndef _CandidateRandomVariables_h
#define _CandidateRandomVariables_h

#include <map>

#include "../RandomVariables/ConstantRandomVariable.h"
#include "../RandomVariables/GaussianRandomVariable.h"
#include "../RandomVariables/AbstractRandomVariable.h"
#include "../Models/AbstractModel.h"

class CandidateRandomVariables {

public:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Key : Name of the parameter. Value : Value of the parameter
    typedef std::map<std::string, double > RandomVariableParameters;              // TODO : Maybe changed?

    /// Key : Name of the random variable. Value : <Type of the candidate random variable, list of parameters>
    typedef std::map<std::string, RandomVariableParameters > RandomVariableParametersMap;

    typedef std::map< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariableMap;


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    CandidateRandomVariables(); // TODO : Add one argument : the XML file
    ~CandidateRandomVariables();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Get a candidate random variable
    std::shared_ptr<AbstractRandomVariable> GetRandomVariable(std::string NameRandomVariable, double Realization);

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the candidate random variables
    void InitializeCandidateRandomVariables(std::shared_ptr<AbstractModel>& M);



protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Gaussian random variable with determined mean
    std::shared_ptr<AbstractRandomVariable> GetGaussianRandomVariable(double Mean, RandomVariableParameters Parameters);

    /// Constant random variable with determined mean
    std::shared_ptr<AbstractRandomVariable> GetConstantRandomVariable(double Mean, RandomVariableParameters Parameters);

    /// Read the parameters of a given random variable within a certain file
    RandomVariableParameters ReadParameters(std::string NameRandomVariable);


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    RandomVariableParametersMap m_RandomVariableParameters;

};


#endif //_CandidateRandomVariables_h
