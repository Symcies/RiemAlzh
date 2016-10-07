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

std::map<std::string, std::vector<double>>
HMWithinGibbsSampler
::Sample(const std::shared_ptr<Realizations>& R, std::shared_ptr<AbstractModel>& M,
         std::shared_ptr<CandidateRandomVariables>& Candidates, const std::shared_ptr<Data>& D)
{
    std::random_device RD;
    std::default_random_engine Generator(RD());
    std::uniform_real_distribution<double> Distribution(0,1);

    std::shared_ptr<Realizations> GibbsRealizations = std::make_shared<Realizations>(*R);
    std::cout << "Realization: ";
    for(Realizations::iterator  it = GibbsRealizations->begin(); it != GibbsRealizations->end(); ++it)
    {
        std::string NameCurrentRV = it->first;
        std::cout << NameCurrentRV << ": ";
        auto CurrentRV = M->GetRandomVariable(NameCurrentRV);

        int i = 0;
        for(std::vector<double>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2, ++i)
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
            
            TestAssert::WarningEquality_Object(0.0, 0.0, "Tau wrong");

            if(log(UnifSample) > Tau) /// It means that the new state is the previous one : no change
            {
                *it2 = CurrentRealization;
                M->UpdateParameters(GibbsRealizations, NameCurrentRV);
                std::cout << CurrentRealization << ". " ;
            }
            else
            {
                //No need to rewrite *it2 as it is already the Candidate Realization
                // *it2 = CandidateRealization;
                std::cout << CandidateRealization << ". ";
            }


            //////////////////////////
            /// DEBUGGING METHODS ////
            //////////////////////////
            /*
            if( NameCurrentRV == "Ksi" or NameCurrentRV == "Tau")
            {
                std::cout << std::endl;
                std::cout << NameCurrentRV << " : " << CurrentRealization << " -> " << CandidateRealization << std::endl;
                std::cout << "Ratio : " << Tau << std::endl;
                std::cout << "Candidate Likelihood/Prior : " << CandidateLogLikelihood << "/" << CandidateLogPrior << std::endl;
                std::cout << "Current   Likelihood/Prior : " << CurrentLogLikelihood << "/" << CurrentLogPrior << std::endl;  
            }
             */
             
            
            /*
            if(true)
            {
                std::cout << NameCurrentRV << " isNaN" << std::endl;
                std::cout << "Candidate Likelihood/Prior : " << CandidateLogLikelihood << "/" << CandidateLogPrior << std::endl;
                std::cout << "Current   Likelihood/Prior : " << CurrentLogLikelihood << "/" << CurrentLogPrior << std::endl;
            }
             */
             
            /*
            if(Tau <10e-4 or Tau > 20)
            {
                std::cout << std::endl << "Ratio of " << NameCurrentRV << " is too small or too large : " << Tau <<  std::endl;
                std::cout << "Candidate Likelihood/Prior : " << CandidateLogLikelihood << "/" << CandidateLogPrior << std::endl;
                std::cout << "Current   Likelihood/Prior : " << CurrentLogLikelihood << "/" << CurrentLogPrior << std::endl;

            }
             */
             
            

            //////////////////////////
            ///   END DEBUGGING   ////
            //////////////////////////

        }
    }
    std::cout << std::endl;

    return *GibbsRealizations;
}
