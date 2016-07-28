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
::Sample(RandomVariable& CurrentRV, double &CurrentRealization, RandomVariable& CandidateRV, AbstractModel &M, Realizations& R, Data& D)
{
    //// Compute the current part
    double CurrentPrior = CurrentRV.second->Likelihood(CurrentRealization);
    double CurrentLikelihood = M.ComputeLikelihood(R, D);
    double Denominator = CurrentPrior * CurrentLikelihood;


    /// Compute the candidate part
    double CandidateRealization = CandidateRV.second->Sample();
    R.at(CurrentRV.first) = CandidateRealization;

    double CandidatePrior = CandidateRV.second->Likelihood(CandidateRealization);
    double CandidateLikelihood = M.ComputeLikelihood(R, D);
    double Numerator = CandidatePrior * CandidateLikelihood;


    // Calculate the hasting metropolis ratio
    double Tau = Numerator/Denominator;

    /// Hasting Metropolis Ratio
    std::random_device RD;
    std::default_random_engine Generator(RD());
    std::uniform_real_distribution<double> Distribution(0,1);
    double unif = Distribution(Generator);

    if(unif > Tau) //The new state is the previous one, no change
    {
        R.at(CurrentRV.first) = CandidateRealization;
    }
}
