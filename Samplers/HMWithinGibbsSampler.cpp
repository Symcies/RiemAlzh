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
::Sample(std::shared_ptr<Realizations>& R, std::shared_ptr<AbstractModel>& M,
         std::shared_ptr<CandidateRandomVariables>& Candidates, const std::shared_ptr<Data>& D)
{
    std::random_device RD;
    std::default_random_engine Generator(RD());
    std::uniform_real_distribution<double> Distribution(0,1);

    for(Realizations::iterator  it = R->begin(); it != R->end(); ++it)
    {
        std::string NameCurrentRV = it->first;
        std::cout << NameCurrentRV << ": ";
        RandomVariable CurrentRV = M->GetRandomVariable(NameCurrentRV);

        int i = 0;
        Realizations R1, R2;
        for(std::vector<double>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2, ++i)
        {
            R1 = *R;

            /// Compute the current part
            std::pair<std::string, double> Realization(NameCurrentRV, i);
            double CurrentRealization = *it2;
            double CurrentPrior = CurrentRV.second->Likelihood(CurrentRealization);
            double CurrentLikelihood = M->ComputeLikelihood(R, D, Realization);

            /// Compute the candidate part
            auto CandidateRV = Candidates->GetRandomVariable(NameCurrentRV, CurrentRealization);
            double CandidateRealization = CandidateRV->Sample();
            double CandidatePrior = CandidateRV->Likelihood(CandidateRealization);
            *it2 = CandidateRealization;
            double CandidateLikelihood = M->ComputeLikelihood(R, D, Realization);

            /// Sampling
            double Tau = CandidatePrior * CandidateLikelihood / (CurrentPrior * CurrentLikelihood);
            // Tau = min(1.0, Tau)
            double UnifSample = Distribution(Generator);

            //////////////////////////
            /// DEBUGGING METHODS ////
            //////////////////////////

            if(std::isnan(Tau))
            {
                std::cout << NameCurrentRV << " isNaN" << std::endl;
                std::cout << "Candidate Likelihood/Prior : " << CandidateLikelihood << "/" << CandidatePrior << std::endl;
                std::cout << "Current   Likelihood/Prior : " << CurrentLikelihood << "/" << CurrentPrior << std::endl;
            }

            if(Tau <10e-4 or Tau > 100)
            {
                std::cout << std::endl << "Ratio of " << NameCurrentRV << " is too small or too large : " << Tau <<  std::endl;
                std::cout << "Candidate Likelihood/Prior : " << CandidateLikelihood << "/" << CandidatePrior << std::endl;
                std::cout << "Current   Likelihood/Prior : " << CurrentLikelihood << "/" << CurrentPrior << std::endl;

            }

            //////////////////////////
            ///   END DEBUGGING   ////
            //////////////////////////

            if(UnifSample > Tau) /// It means that the new state is the previous one : no change
            {
                *it2 = CurrentRealization;
                std::cout << "0. " ;
            }
            else
            {
                // TO DO : Should be in the compute likelihood part ?
                // Useless if there are only pointers to the realizations
                std::cout << "1. ";
            }

            R2 = *R;
            if(R1 == R2)
            {
                //std::cout << "Donc ca marche? => " << UnifSample <<  ">" <<  Tau << std::endl;
            }
        }
    }
    std::cout << std::endl;
}
