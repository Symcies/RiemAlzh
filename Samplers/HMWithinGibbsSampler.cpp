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
        RandomVariable CurrentRV = M->GetRandomVariable(NameCurrentRV);

        for(std::vector<double>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {

            /// Compute the current part
            double CurrentRealization = *it2;
            double CurrentPrior = CurrentRV.second->Likelihood(CurrentRealization);
            double CurrentLikelihood = M->ComputeLikelihood(R, D);

            /// Compute the candidate part
            std::shared_ptr<AbstractRandomVariable> CandidateRV = M->GetCandidateRandomVariable(NameCurrentRV, CurrentRealization);
            double CandidateRealization = CandidateRV->Sample();
            double CandidatePrior = CandidateRV->Likelihood(CandidateRealization);
            // TODO : gros problème car ça omet les realizations deltak qui ne sont updatées qu'avec setrealization
            *it2 = CandidateRealization;
            double CandidateLikelihood = M->ComputeLikelihood(R, D);

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
                std::cout << "Likelihood : " << CandidateLikelihood << std::endl;
                //std::cout << "Likelihood : " << CurrentLikelihood << std::endl;
            }

            //////////////////////////
            ///   END DEBUGGING   ////
            //////////////////////////

            std::cout << Tau << ". ";
            if(UnifSample > Tau) /// It means that the new state is the previous one : no change
            {
                *it2 = CurrentRealization;
                //std::cout << "0. " ;
            }
            else
            {
                // TO DO : Should be in the compute likelihood part ?
                // Useless if there are only pointers to the realizations
                M->SetRealization(NameCurrentRV, CandidateRealization);
                //std::cout << "1. ";
            }


        }
    }
    std::cout << std::endl;
}
