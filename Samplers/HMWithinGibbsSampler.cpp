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
    m_CandidateRandomVariables.InitializeCandidateRandomVariables(R);
}


HMWithinGibbsSampler::MultiRealizations
HMWithinGibbsSampler
::Sample(const std::shared_ptr<MultiRealizations>& R, std::shared_ptr<AbstractModel>& M,
         const std::shared_ptr<Data>& D, int IterationNumber)
{
    std::random_device RD;
    std::mt19937 Generator(RD());
    std::uniform_real_distribution<double> Distribution(0.0, 1.0);

    std::shared_ptr<MultiRealizations> GibbsRealizations = std::make_shared<MultiRealizations>(*R);
    for(auto&& it : *GibbsRealizations)
    {
        std::string NameCurrentRV = it.first;
        auto CurrentRV = M->GetRandomVariable(NameCurrentRV);

        int i = 0;
        for(auto it2 = it.second.begin(); it2 != it.second.end(); ++it2, ++i)
        {

            /// Compute the current part
            double CurrentRealization = *it2;
            double CurrentLogPrior = CurrentRV->LogLikelihood(CurrentRealization);
            double CurrentLogLikelihood = M->ComputeLogLikelihood(GibbsRealizations, D, std::pair<std::string, int> (NameCurrentRV, i));

            /// Compute the candidate part
            auto CandidateRV = m_CandidateRandomVariables.GetRandomVariable(NameCurrentRV, i, CurrentRealization);
            double CandidateRealization = CandidateRV.Sample();
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
            }
            else
            {
                //No need to rewrite *it2 as it is already the Candidate Realization
            }
            
            //UpdatePropositionDistributionVariance(CandidateRV, exp(Tau), IterationNumber);
            
            
        }
    }

    return *GibbsRealizations;
}
