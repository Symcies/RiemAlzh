#include "HMWithinGibbsSampler.h"

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
::Sample(std::shared_ptr<Realizations>& R, std::shared_ptr<AbstractModel>& M, const std::shared_ptr<Data>& D)
{
    std::random_device RD;
    std::default_random_engine Generator(RD());
    std::uniform_real_distribution<double> Distribution(0,1);

    for(Realizations::iterator  it = R->begin(); it != R->end(); ++it)
    {
        std::string NameCurrentRV = it->first;
        std::cout << std::endl << NameCurrentRV << std::endl;
        RandomVariable CurrentRV = M->GetRandomVariable(NameCurrentRV);

        for(std::vector<double>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
            /// Compute the current part
            double CurrentRealization = *it2;
            double CurrentPrior = CurrentRV.second->Likelihood(CurrentRealization);
            double CurrentLikelihood = M->ComputeLikelihood(R, D);

            /// Compute the candidate part
            ////// TODO : GET THE CANDIDATE PART AND CHANGE THE FOLLOWING LINES
            GaussianRandomVariable CandidateRV = GaussianRandomVariable(CurrentRealization, 0.01);
            double CandidateRealization = CandidateRV.Sample();
            //////////////////////////////
            double CandidatePrior = CandidateRV.Likelihood(CandidateRealization);
            *it2 = CandidateRealization;
            double CandidateLikelihood = M->ComputeLikelihood(R, D);

            /// Sampling
            double Tau = CandidatePrior * CandidateLikelihood / (CurrentPrior * CurrentLikelihood);
            double UnifSample = Distribution(Generator);

            if(UnifSample > Tau) /// It means that the new state is the previous one : no change
            {
                *it2 = CurrentRealization;
            }
            else
            {
                M->SetRealization(NameCurrentRV, CandidateRealization);
            }


        }
    }
}
