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
        auto CurrentRV = M->GetRandomVariable(NameCurrentRV);

        int i = 0;
        for(std::vector<double>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2, ++i)
        {

            /// Compute the current part
            double CurrentRealization = *it2;
            double CurrentPrior = CurrentRV->Likelihood(CurrentRealization);
            double CurrentLogLikelihood = M->ComputeLogLikelihood(R, D, std::pair<std::string, int> (NameCurrentRV, i));

            /// Compute the candidate part
            auto CandidateRV = Candidates->GetRandomVariable(NameCurrentRV, CurrentRealization);
            double CandidateRealization = CandidateRV->Sample();
            double CandidatePrior = CurrentRV->Likelihood(CandidateRealization);
            *it2 = CandidateRealization;
            double CandidateLogLikelihood = M->ComputeLogLikelihood(R, D, std::pair<std::string, int> (NameCurrentRV, i));

            /// Sampling
            double Tau = CandidatePrior * CandidateLogLikelihood / (CurrentPrior * CurrentLogLikelihood);
            double UnifSample = Distribution(Generator);


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


            //////////////////////////
            /// DEBUGGING METHODS ////
            //////////////////////////

            if(std::isnan(Tau))
            {
                std::cout << NameCurrentRV << " isNaN" << std::endl;
                std::cout << "Candidate Likelihood/Prior : " << CandidateLogLikelihood << "/" << CandidatePrior << std::endl;
                std::cout << "Current   Likelihood/Prior : " << CurrentLogLikelihood << "/" << CurrentPrior << std::endl;
            }

            if(Tau <10e-4 or Tau > 20)
            {
                std::cout << std::endl << "Ratio of " << NameCurrentRV << " is too small or too large : " << Tau <<  std::endl;
                std::cout << "Candidate Likelihood/Prior : " << CandidateLogLikelihood << "/" << CandidatePrior << std::endl;
                std::cout << "Current   Likelihood/Prior : " << CurrentLogLikelihood << "/" << CurrentPrior << std::endl;

            }

            //////////////////////////
            ///   END DEBUGGING   ////
            //////////////////////////

        }
    }
    std::cout << std::endl;
}
