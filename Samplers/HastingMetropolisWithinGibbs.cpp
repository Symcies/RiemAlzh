#include "HastingMetropolisWithinGibbs.h"


HastingMetropolisWithinGibbs
::HastingMetropolisWithinGibbs()
{

}


bool
HastingMetropolisWithinGibbs
::Sample(AbstractRandomVariable &CurrentRandomVariable, AbstractRandomVariable &CandidateRandomVariable, LongitudinalModel &LM)
{


    double CurrentState = CurrentRandomVariable.GetCurrentState();

    // Calculate the law for the current state
    double CurrentDensity = CurrentRandomVariable.GetDensity();
    double CurrentLikelihood = LM.GetCurrentLikelihood();
    double RatioDenominator = CurrentDensity* CurrentLikelihood;

    std::cout << " Current Density / Likelihood : " << CurrentDensity << " / " << CurrentLikelihood << std::endl;


    // Change the manifold
    CurrentRandomVariable.SetCurrentState(CandidateRandomVariable.GetCurrentState());

    // Calculate the law of the candidate state
    double CandidateDensity = CandidateRandomVariable.GetDensity();
    double CandidateLikelihood = LM.ComputeLikelihood();
    double RatioNumerator = CurrentRandomVariable.GetDensity()* CandidateLikelihood;

    std::cout << " Candidate Density / Likelihood : " << CandidateDensity << " / " << CandidateLikelihood << std::endl;

    // Calculate the hasting metropolis ratio
    double Tau = RatioNumerator/RatioDenominator;


    //std::cout << "Current Density: " << CurrentDensity << " - Current Likelihood : " << LM.GetCurrentLikelihood() << " - Numerator Ratio : " << RatioNumerator << std::endl;
    //std::cout << "Candida Density: " << CurrentRandomVariable.GetDensity() << " - Candida Likelihood : " << CandidateLikelihood << " - Denom Ratio : " << RatioDenominator << std::endl;

    /// Hasting Metropolis Ratio
    std::random_device RD;
    std::default_random_engine Generator(RD());
    std::uniform_real_distribution<double> Distribution(0,1);
    double unif = Distribution(Generator);

    // If the new state is rejected, then the manifold has to be set to its initial form
    if(unif > Tau)    // The new state is the previous one : no change
    {
        CurrentRandomVariable.SetCurrentState(CurrentState);
        return false;
    }
    else             // The new state is the candidate state
    {
        LM.SetCurrentLikelihood(CandidateLikelihood);
        return true;
    }


}
