#include "BlockedGibbsSampler2.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

BlockedGibbsSampler2
::BlockedGibbsSampler2() 
{
    
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void
BlockedGibbsSampler2
::InitializeSampler(const std::shared_ptr<MultiRealizations> &R) 
{
    m_CandidateRandomVariables.InitializeCandidateRandomVariables(R);
    
    for(auto it = R->begin(); it != R->end(); ++it)
    {
        if(it->second.size() == 1)
        {
            m_NamePopulationVariables.push_back(it->first);
            std::cout << "Population: " << it->first << std::endl;
        }
        else
        {
            m_NameIndividualVariables.push_back(it->first);
            std::cout << "Individual: " << it->first << std::endl;
        }
    }
}

BlockedGibbsSampler2::MultiRealizations
BlockedGibbsSampler2
::Sample(const std::shared_ptr<MultiRealizations> &R, std::shared_ptr<AbstractModel> &M,
         const std::shared_ptr<Data> &D, int IterationNumber) 
{

    /// Sample the population-wide random variables
    auto NewRealizations = std::make_shared<MultiRealizations>( SamplePopulation(R, M, D) );
    
    /// Sample the candidate random variables
    for(int i = 0; i < D->size(); ++i)
    {
        *NewRealizations = SampleIndividual(i, NewRealizations, M, D);
    }
    
    // TODO : Add the adaptive step within the two samplings : need to be done when the ratio is calculated
    
    return *NewRealizations;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

BlockedGibbsSampler2::MultiRealizations
BlockedGibbsSampler2
::SamplePopulation(const std::shared_ptr<MultiRealizations> &R, std::shared_ptr<AbstractModel> &M, const std::shared_ptr<Data> &D)
{
    /// Get Z_pop(k) : the current population realizations
    auto CurrentPopRealizations = std::make_shared<UniqueRealizations> (GetCurrentPopulationRealizations(R) );
    
    /// Get Z(*) : the candidate population realizations
    UniqueRealizations CandidatePopRealizations = GetCandidateRealizations(CurrentPopRealizations, 0);
    
    /// Get the acceptance ratio
    double Ratio = - M->ComputeLogLikelihood(R, D);
    for(const auto& it : *CurrentPopRealizations)
    {
        auto RandomVariable = M->GetRandomVariable(it.first);
        Ratio -= RandomVariable->LogLikelihood(it.second);
        double CandidateRealization = CandidatePopRealizations.at(it.first);
        Ratio += RandomVariable->LogLikelihood( CandidateRealization );
    }
    
    auto NewRealization = std::make_shared<MultiRealizations>(*R);
    for(const auto& it : CandidatePopRealizations)
    {
        NewRealization->at(it.first) = VectorType(1, it.second );
    }
    Ratio += M->ComputeLogLikelihood(NewRealization, D);
    
    /// Compare the ratio to a uniform distribution
    std::random_device RD;
    std::mt19937 Generator(RD());
    std::uniform_real_distribution<double> Distribution(0.0, 1.0);
    double UniformSample = Distribution(Generator);
    
    if(log(UniformSample) > Ratio)
    {
        return *R;
    }
    else
    {
        return *NewRealization;
    }
}


std::map<std::string, double>
BlockedGibbsSampler2
::GetCurrentPopulationRealizations(const std::shared_ptr<MultiRealizations> &R) 
{
    UniqueRealizations PopReal;
    for(const auto& it : m_NamePopulationVariables)
    {
        PopReal[it] = R->at(it)(0);
    }
    
    return PopReal;
}

std::map<std::string, double>
BlockedGibbsSampler2
::GetCurrentIndividualRealizations(const std::shared_ptr<MultiRealizations> &R, int i) 
{
    UniqueRealizations IndivReal;
    for(const auto& it : m_NameIndividualVariables)
    {
        IndivReal[it] = R->at(it)(i);
    }
    return IndivReal;
}

std::map<std::string, double>
BlockedGibbsSampler2
::GetCandidateRealizations(const std::shared_ptr<UniqueRealizations> &R, int Number) 
{
    // The Number argument corresponds to the realization number of the random variable
    // If Population Random Variable, then Number = 0
    // If Individual Random Variables, then Number belongs to [0, NumberOfSubjects]
    
    
    UniqueRealizations CandidateRealizations;
    for(const auto& it : *R)
    {
        double CurrentRealization = it.second;
        auto CandidateRandomVariable = m_CandidateRandomVariables.GetRandomVariable( it.first, Number, CurrentRealization );
        double CandidateRealization = CandidateRandomVariable.Sample();
        
        CandidateRealizations[it.first] = CandidateRealization;
    }

    return CandidateRealizations;
}

BlockedGibbsSampler2::MultiRealizations
BlockedGibbsSampler2
::SampleIndividual(int i, const std::shared_ptr<MultiRealizations> &R, std::shared_ptr<AbstractModel> &M,
                   const std::shared_ptr<Data> &D) 
{
    /// Get the current realization Z(k)_i of individual i
    auto CurrentIndivRealizations = std::make_shared<UniqueRealizations>( GetCurrentIndividualRealizations(R, i) ); 
    
    /// Get the candidate realization Z*_i of the individual i
    UniqueRealizations CandidateIndivRealizations = GetCandidateRealizations(CurrentIndivRealizations, i);
    
    ///Get the acceptance ratio
    double Ratio = - M->ComputeLogLikelihood(R, D);
    for(const auto& it : *CurrentIndivRealizations)
    {
        auto RandomVariable = M->GetRandomVariable(it.first);
        Ratio -= RandomVariable->LogLikelihood(it.second);
        double CandidateRealization = CandidateIndivRealizations.at(it.first);
        Ratio += RandomVariable->LogLikelihood( CandidateRealization );
    }
    
    auto NewRealization = std::make_shared<MultiRealizations>(*R);
    for(const auto& it : CandidateIndivRealizations)
    {
        NewRealization->at(it.first)(i) = it.second;
    }
    Ratio += M->ComputeLogLikelihood(NewRealization, D);
    
    /// Compare the ratio to a uniform distribution
    std::random_device RD;
    std::mt19937 Generator(RD());
    std::uniform_real_distribution<double> Distribution(0.0, 1.0);
    double UniformSample = Distribution(Generator);
    
    if(log(UniformSample) > Ratio)
    {
        return *R;
    }
    else
    {
        return *NewRealization;
    }
}