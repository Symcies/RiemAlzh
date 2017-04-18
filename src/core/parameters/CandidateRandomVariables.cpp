#include <stdexcept>
#include "CandidateRandomVariables.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


CandidateRandomVariables::CandidateRandomVariables()
{

}

CandidateRandomVariables::~CandidateRandomVariables()
{
    // TODO : How to generalize?
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

GaussianRandomVariable CandidateRandomVariables::GetRandomVariable(std::string rand_var_name,int num_real)
const
{
    return new_proposition_distribution_.at(string_to_int_key_.at(rand_var_name)).at(num_real);
}


GaussianRandomVariable CandidateRandomVariables::GetRandomVariable(int rand_var_key, int num_real)
const
{
    return new_proposition_distribution_.at(rand_var_key).at(num_real);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void CandidateRandomVariables::InitializeCandidateRandomVariables(const Realizations& reals, const AbstractModel& model)
{

    for(auto it = reals.begin(); it != reals.end(); ++it)
    {
        std::vector< GaussianRandomVariable > prop_distrib_per_real;
        std::string name = reals.ReverseKeyToName(it->first);
        int_to_string_key_.insert({it->first, name});
        string_to_int_key_.insert({name, it->first});

        for(auto it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
            double variance = model.InitializePropositionDistributionVariance(name);
            GaussianRandomVariable gaussian_rand_var(*it2, variance);
            prop_distrib_per_real.push_back( gaussian_rand_var );
        }

        new_proposition_distribution_[it->first] = prop_distrib_per_real;
    }

}


void CandidateRandomVariables::UpdatePropositionVariableVariance(std::string name, int num_real, ScalarType new_variance)
{

    new_proposition_distribution_.at(string_to_int_key_.at(name))[num_real].Update({{"Variance", new_variance}});
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////
