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
    auto GibbsRealizations = std::make_shared<Realizations>(*R);
    double CurrentLogLikelihood = M->ComputeLogLikelihood(R, D);
    
    /// Loop over each realization
    for(auto it : *R)
    {
        /// Get the random variable associated to the realizations
        std::string NameCurrentVariable = it.first;
        auto CurrentRandomVariable = M->GetRandomVariable( NameCurrentVariable );
        
        std::cout << NameCurrentVariable << " : ";
        /// Initialize loop
        auto IterReal = it.second.begin();
        int i = 0;
        for( ; IterReal != it.second.end(); ++IterReal, ++i)
        {
            
            /// Current random vabiable 
            double CurrentRealization = *IterReal;
            double CurrentLogPrior = CurrentRandomVariable->LogLikelihood(*IterReal);
            
            /// Candidate random variable
            auto CandidateRandomVariable = Candidates->GetRandomVariable(NameCurrentVariable, *IterReal);
            double CandidateRealization = CandidateRandomVariable->Sample();
            double CandidateLogPrior = CurrentRandomVariable->LogLikelihood( CandidateRealization );
            *IterReal = CandidateRealization;
            double CandidateLogLikelihood = M->ComputeLogLikelihood(R, D);
            
            /// Sampling
            double Tau = CandidateLogPrior + CandidateLogLikelihood - (CurrentLogPrior + CurrentLogLikelihood);
            Tau = std::min(0.0, Tau);
            double UniformSampler = Distribution(Generator);
            
            /// Acceptation / Reject
            if(log(UniformSampler) > Tau)   /// Reject : no changes
            {
                std::cout << CurrentRealization << ". " ;
            }
            else                            /// Acceptation
            {
                GibbsRealizations->at(NameCurrentVariable)[i] = *IterReal;
                std::cout << CandidateRealization << ". ";
            }
            
            *IterReal = CurrentRealization;
        }
    }
    std::cout << std::endl;
    
    return *GibbsRealizations;
}