#include <stdexcept>
#include "CandidateRandomVariables.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


CandidateRandomVariables
::CandidateRandomVariables()
{

}

CandidateRandomVariables
::~CandidateRandomVariables()
{
    // TODO : How to generalize?
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<AbstractRandomVariable>
CandidateRandomVariables
::GetRandomVariable(std::string NameRandomVariable, double Realization)
{
    if(NameRandomVariable.size() != 1 and NameRandomVariable.find_first_of("#") != std::string::npos)
    {
        NameRandomVariable = NameRandomVariable.substr(0, NameRandomVariable.find_first_of("#"));
    }
    

    RandomVariableParameters Parameters = m_RandomVariableParameters.at(NameRandomVariable);

    double RVType = Parameters.at("Type");

    if(RVType == 1)
    {
        return GetGaussianRandomVariable(Realization, Parameters);
    }
    else if(RVType == 2)
    {
        return GetConstantRandomVariable(Realization, Parameters);
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
CandidateRandomVariables
::InitializeCandidateRandomVariables(std::shared_ptr<AbstractModel> &M)
{
    /// Initialization
    RandomVariableParametersMap MapParameters;


    RandomVariableMap ModelRandomVariables = M->GetRandomVariables();
    for(auto it = ModelRandomVariables.begin(); it != ModelRandomVariables.end(); ++it)
    {
        std::string NameRandomVariable = it->first;
        /// It means that a parameter is shared among different random variables
        // TODO : Write somewhere what the # is for : to explain that there is something shared
        // For instance : Delta0, Delta 1, ... => Delta#1, Delta#2 , ... because sharing the variance!
        if(NameRandomVariable.size() != 1 and NameRandomVariable.find_first_of("#") != std::string::npos)
        {
            NameRandomVariable = NameRandomVariable.substr(0, NameRandomVariable.find_first_of("#"));
        }

        MapParameters[NameRandomVariable] = ReadParameters(NameRandomVariable);
    }

    m_RandomVariableParameters = MapParameters;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<AbstractRandomVariable>
CandidateRandomVariables
::GetGaussianRandomVariable(double Mean, RandomVariableParameters Parameters)
{
    double Variance = Parameters.at("Variance");
    std::shared_ptr<AbstractRandomVariable> RandomVariable = std::make_shared<GaussianRandomVariable>(Mean, Variance);
    return RandomVariable;
}


std::shared_ptr<AbstractRandomVariable>
CandidateRandomVariables
::GetConstantRandomVariable(double Mean, RandomVariableParameters Parameters) 
{
    std::shared_ptr<AbstractRandomVariable> RandomVariable = std::make_shared<ConstantRandomVariable>(Mean);
    return RandomVariable;
}



std::map<std::string, double >
CandidateRandomVariables
::ReadParameters(std::string NameRandomVariable)
{


    /////////////////////////////////////////////////////
    /// The following has to follow functions
    /// e.g. : Types = READ(XMLFile)
    /// e.g. : Parameters = READ(XMLFile)
    /////////////////////////////////////////////////////

    /// QUICK CHANGES IN THE INITIALIZATION
    double P0Variance = 0.00000005;
    double T0Variance = 0.000005;
    double V0Variance = 0.000000001;
    double DeltaVariance = 0.000005;
    double BetaVariance = 0.00000002;
    double TauVariance = 0.8;
    double KsiVariance = 0.003;
    double SVariance = 0.5;



    RandomVariableParameters Parameters;

    if(NameRandomVariable == "P0")
    {
        Parameters["Type"] = 1;
        Parameters["Variance"] = P0Variance;
    }
    else if(NameRandomVariable == "T0")
    {
        Parameters["Type"] = 1;
        Parameters["Variance"] = T0Variance;
    }
    else if(NameRandomVariable == "V0")
    {
        Parameters["Type"] = 1;
        Parameters["Variance"] = V0Variance;
    }
    else if(NameRandomVariable == "Tau")
    {
        Parameters["Type"] = 1;
        Parameters["Variance"] = TauVariance;
    }
    else if(NameRandomVariable == "Ksi")
    {
        Parameters["Type"] = 1;
        Parameters["Variance"] = KsiVariance;
    }
    else if(NameRandomVariable == "Delta")
    {
        Parameters["Type"] = 1;
        Parameters["Variance"] = DeltaVariance;
    }
    else if(NameRandomVariable == "Beta")
    {
        Parameters["Type"] = 1;
        Parameters["Variance"] = BetaVariance;
    }
    else if(NameRandomVariable == "S")
    {
        Parameters["Type"] = 1;
        Parameters["Variance"] = SVariance;
    }
    else
    {
        std::cout << NameRandomVariable << std::endl;
        throw std::invalid_argument("This random variables - defined in the model - does not exist in the XML File");
    }


    /////////////////////////////////////////////////////
    /// End of the future XMLReader function
    /////////////////////////////////////////////////////

    return Parameters;
}