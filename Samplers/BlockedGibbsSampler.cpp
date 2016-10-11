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

std::map<std::string, std::vector<double>>
BlockedGibbsSampler
::Sample(const std::shared_ptr<Realizations> &R, std::shared_ptr<AbstractModel> &M,
         std::shared_ptr<CandidateRandomVariables> &Candidates, const std::shared_ptr<Data> &D) 
{
    
    /// Initialize Random number generator for the acceptation / rejection part
    std::random_device RD;
    std::default_random_engine Generator(RD());
    std::uniform_real_distribution<double> Distribution(0, 1);
    
    /// Initialize the realization to return and the corresponding likelihood
    
    auto OutputRealizations = std::make_shared<Realizations>(*R);
    auto WorkingRealizations = std::make_shared<Realizations>(*R);
    double CurrentLogLikelihood = M->ComputeLogLikelihood(R, D);
    
    /// Loop over each realization
    for(auto it = WorkingRealizations->begin(); it != WorkingRealizations->end(); ++it)
    {
        /// Get the random variable associated to the realizations
        std::string NameCurrentVariable = it->first;
        auto CurrentRandomVariable = M->GetRandomVariable( NameCurrentVariable );
        
        std::cout << NameCurrentVariable << " : ";
        /// Initialize loop
        
        int i = 0;
        for(auto IterReal = it->second.begin() ; IterReal != it->second.end(); ++IterReal, ++i)
        {
            
            /// Current random vabiable 
            double CurrentRealization = *IterReal;
            double CurrentLogPrior = CurrentRandomVariable->LogLikelihood(*IterReal);
            
            /// Candidate random variable
            auto CandidateRandomVariable = Candidates->GetRandomVariable(NameCurrentVariable, CurrentRealization );
            double CandidateRealization = CandidateRandomVariable->Sample();
            double CandidateLogPrior = CurrentRandomVariable->LogLikelihood( CandidateRealization );
            *IterReal = CandidateRealization;
            double CandidateLogLikelihood = M->ComputeLogLikelihood(WorkingRealizations, D);
            
            /// Sampling
            double Tau = CandidateLogPrior + CandidateLogLikelihood - (CurrentLogPrior + CurrentLogLikelihood);
            Tau = std::min(0.0, Tau);
            double UniformSampler = Distribution(Generator);
            
            /// Acceptation / Reject
            if(log(UniformSampler) > Tau)   /// Reject : no changes
            {
                OutputRealizations->at(NameCurrentVariable)[i] = CurrentRealization;
                std::cout << CurrentRealization << ". " ;
            }
            else                            /// Acceptation
            {
                OutputRealizations->at(NameCurrentVariable)[i] = CandidateRealization;
                std::cout << CandidateRealization << ". ";
            }
            
            *IterReal = CurrentRealization;
            M->UpdateParameters(R, NameCurrentVariable);
            
        }
    }
    std::cout << std::endl;
    
    return *OutputRealizations;
}