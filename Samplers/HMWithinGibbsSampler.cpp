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
::Sample(RandomVariable& CurrentRV, std::shared_ptr< AbstractRandomVariable >& CandidateRV, std::shared_ptr<AbstractModel>& M, Realizations& R, const std::shared_ptr<Data>& D)
{
    //// Compute the current part
    double CurrentRealization = R.at(CurrentRV.first);
    double CurrentPrior = CurrentRV.second->Likelihood(CurrentRealization);
    double CurrentLikelihood = M->ComputeLikelihood(R, D);
    double Denominator = CurrentPrior * CurrentLikelihood;


    /// Compute the candidate part
    double CandidateRealization = CandidateRV->Sample();
    R.at(CurrentRV.first) = CandidateRealization;

    double CandidatePrior = CandidateRV->Likelihood(CandidateRealization);
    double CandidateLikelihood = M->ComputeLikelihood(R, D);
    double Numerator = CandidatePrior * CandidateLikelihood;


    // Calculate the hasting metropolis ratio
    double Tau = Numerator/Denominator;


    /// Hasting Metropolis Ratio
    std::random_device RD;
    std::default_random_engine Generator(RD());
    std::uniform_real_distribution<double> Distribution(0,1);
    double Unif = Distribution(Generator);

    if(Unif > Tau) //The new state is the previous one, no change
    {
        R.at(CurrentRV.first) = CandidateRealization;
    }

    //std::cout << "Name : " << CurrentRV.first << std::endl;
    //std::cout << "Candidate Prior/Likelihood : "<< CandidatePrior << "/" << CandidateLikelihood << std::endl;
    //std::cout << "Current Prior/Likelihood : "<< CurrentPrior << "/" << CurrentLikelihood << std::endl;
}
