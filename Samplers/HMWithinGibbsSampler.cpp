#include "HMWithinGibbsSampler.h"
#include "../Tests/TestAssert.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

HMWithinGibbsSampler
::HMWithinGibbsSampler() 
{

}




////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
HMWithinGibbsSampler
::InitializeSampler(const std::shared_ptr<MultiRealizations> &R) 
{
    
}


std::map<std::string, std::vector<double>>
HMWithinGibbsSampler
::Sample(const std::shared_ptr<MultiRealizations>& R, std::shared_ptr<AbstractModel>& M,
         std::shared_ptr<CandidateRandomVariables>& Candidates, const std::shared_ptr<Data>& D)
{
    std::random_device RD;
    std::mt19937 Generator(RD());
    std::uniform_real_distribution<double> Distribution(0.0, 1.0);

    std::shared_ptr<MultiRealizations> GibbsRealizations = std::make_shared<MultiRealizations>(*R);
    std::cout << "Realization: ";
    for(auto&& it : *GibbsRealizations)
    {
        std::string NameCurrentRV = it.first;
        std::cout << NameCurrentRV << ": ";
        auto CurrentRV = M->GetRandomVariable(NameCurrentRV);

        int i = 0;
        for(auto it2 = it.second.begin(); it2 != it.second.end(); ++it2, ++i)
        {

            /// Compute the current part
            double CurrentRealization = *it2;
            double CurrentLogPrior = CurrentRV->LogLikelihood(CurrentRealization);
            double CurrentLogLikelihood = M->ComputeLogLikelihood(GibbsRealizations, D, std::pair<std::string, int> (NameCurrentRV, i));

            /// Compute the candidate part
            auto CandidateRV = Candidates->GetRandomVariable(NameCurrentRV, CurrentRealization);
            double CandidateRealization = CandidateRV->Sample();
            double CandidateLogPrior = CurrentRV->LogLikelihood(CandidateRealization);
            *it2 = CandidateRealization;
            double CandidateLogLikelihood = M->ComputeLogLikelihood(GibbsRealizations, D, std::pair<std::string, int> (NameCurrentRV, i));
            
            /// Sampling
            double Tau = CandidateLogPrior + CandidateLogLikelihood - (CurrentLogPrior + CurrentLogLikelihood);
            Tau = std::min(0.0, Tau);
            double UnifSample = Distribution(Generator);

            if(log(UnifSample) > Tau) /// It means that the new state is the previous one : no change
            {
                *it2 = CurrentRealization;
                M->UpdateParameters(GibbsRealizations, NameCurrentRV);
                std::cout << CurrentRealization << ". " ;
            }
            else
            {
                //No need to rewrite *it2 as it is already the Candidate Realization
                std::cout << CandidateRealization << ". ";
            }
        }
    }
    std::cout << std::endl;

    return *GibbsRealizations;
}
