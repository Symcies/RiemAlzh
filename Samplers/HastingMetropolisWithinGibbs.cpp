#include "HastingMetropolisWithinGibbs.h"


HastingMetropolisWithinGibbs
::HastingMetropolisWithinGibbs()
{

}




void
HastingMetropolisWithinGibbs
::Sample(AbstractRandomVariable &CurrentRandomVariable, AbstractRandomVariable &CandidateRandomVariable, LongitudinalModel &LM)
{


    double CurrentState = CurrentRandomVariable.GetCurrentState();

    // Calculate the law for the current state
    double RatioDenominator = CurrentRandomVariable.GetDensity()*LM.ComputeLikelihood();

    // Draw a candidate state
    double CandidateVariance = 1.0;
    auto Candidate = GaussianRandomVariable(CurrentState, CandidateVariance);

    // Change the manifold
    CurrentRandomVariable.SetCurrentState(Candidate.GetCurrentState());

    // Calculate the law of the candidate state
    double RatioNumerator = CurrentRandomVariable.GetDensity()* LM.ComputeLikelihood();

    // Calculate the hasting metropolis ratio
    double Tau = RatioNumerator/RatioDenominator;

    /// Hasting Metropolis Ratio
    std::random_device RD;
    std::default_random_engine Generator(RD());
    std::uniform_real_distribution<double> Distribution(0,1);
    double unif = Distribution(Generator);

    // If the new state is rejected, then the manifold has to be set to its initial form
    if(unif > Tau)
    {
        CurrentRandomVariable.SetCurrentState(CurrentState);
    }

}
