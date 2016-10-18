#include "BlockedGibbsSampler.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

BlockedGibbsSampler
::BlockedGibbsSampler() 
{
    
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void
BlockedGibbsSampler
::InitializeSampler(const std::shared_ptr<MultiRealizations> &R) 
{
    for(auto it = R->begin(); it != R->end(); ++it)
    {
        if(it->second.size() == 1)
        {
            m_NamePopulationVariables.push_back(it->first);
        }
        else
        {
            m_NameIndividualVariables.push_back(it->first);
        }
    }
}

std::map<std::string, std::vector<double>>
BlockedGibbsSampler
::Sample(const std::shared_ptr<MultiRealizations> &R, std::shared_ptr<AbstractModel> &M,
         std::shared_ptr<CandidateRandomVariables> &Candidates, const std::shared_ptr<Data> &D) 
{

    auto NewRealizations = std::make_shared<MultiRealizations>( SamplePopulation(R, M, Candidates, D) );
    
    for(int i = 0; i < D->size(); ++i)
    {
        *NewRealizations = SampleIndividual(i, NewRealizations, M, Candidates, D);
    }
    
    return *NewRealizations;
    
    
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

std::map<std::string, std::vector<double>>
BlockedGibbsSampler
::SamplePopulation(const std::shared_ptr<MultiRealizations> &R, std::shared_ptr<AbstractModel> &M,
                   std::shared_ptr<CandidateRandomVariables> &Candidates, const std::shared_ptr<Data> &D)
{
    /// Get Z_pop(k) : the current population realizations
    auto CurrentPopRealizations = std::make_shared<UniqueRealizations> (GetCurrentPopulationRealizations(R) );
    
    /// Get Z(*) : the candidate population realizations
    UniqueRealizations CandidatePopRealizations = GetCandidateRealizations(CurrentPopRealizations, Candidates);
    
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
        NewRealization->at(it.first) = { it.second };
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
BlockedGibbsSampler
::GetCurrentPopulationRealizations(const std::shared_ptr<MultiRealizations> &R) 
{
    UniqueRealizations PopReal;
    for(const auto& it : m_NamePopulationVariables)
    {
        PopReal[it] = R->at(it)[0];
    }
    
    return PopReal;
}

std::map<std::string, double>
BlockedGibbsSampler
::GetCurrentIndividualRealizations(const std::shared_ptr<MultiRealizations> &R, int i) 
{
    UniqueRealizations IndivReal;
    for(const auto& it : m_NameIndividualVariables)
    {
        IndivReal[it] = R->at(it)[i];
    }
    return IndivReal;
}

std::map<std::string, double>
BlockedGibbsSampler
::GetCandidateRealizations(const std::shared_ptr<UniqueRealizations> &R,
                                     const std::shared_ptr<CandidateRandomVariables> &Candidates) 
{
    UniqueRealizations CandidateRealizations;
    for(const auto& it : *R)
    {
        double CurrentRealization = R->at(it.first);
        auto CandidateRandomVariable = Candidates->GetRandomVariable( it.first, CurrentRealization );
        double CandidateRealization = CandidateRandomVariable->Sample();
        
        CandidateRealizations[it.first] = CandidateRealization;
    }

    return CandidateRealizations;
}

std::map<std::string, std::vector<double>> 
BlockedGibbsSampler
::SampleIndividual(int i, const std::shared_ptr<MultiRealizations> &R, std::shared_ptr<AbstractModel> &M,
                   std::shared_ptr<CandidateRandomVariables> &Candidates,
                   const std::shared_ptr<Data> &D) 
{
    /// Get the current realization Z(k)_i of individual i
    auto CurrentIndivRealizations = std::make_shared<UniqueRealizations>( GetCurrentIndividualRealizations(R, i) ); 
    
    /// Get the candidate realization Z*_i of the individual i
    UniqueRealizations CandidateIndivRealizations = GetCandidateRealizations(CurrentIndivRealizations, Candidates);
    
    ///Get the acceptance ratio
    double Ratio = - M->ComputeIndividualLogLikelihood(R, D, i);
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
        NewRealization->at(it.first)[i] = it.second;
    }
    Ratio += M->ComputeIndividualLogLikelihood(NewRealization, D, i);
    
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