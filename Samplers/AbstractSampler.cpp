#include "AbstractSampler.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


AbstractSampler
::AbstractSampler()
{

}

AbstractSampler
::~AbstractSampler()
{

}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Methods :
////////////////////////////////////////////////////////////////////////////////////////////////////

double
AbstractSampler
::DecreasingStepSize(int Iteration, int NoMemoryTime) 
{
    double Epsilon = std::max(1, Iteration - NoMemoryTime);
    double Lambda = 0.6; // Lambda should belong to ]1/2 ; 1]
    return 1.0 / pow(Epsilon, Lambda);
}

void
AbstractSampler
::UpdatePropositionDistributionVariance(GaussianRandomVariable &GRV, double Ratio, int Iteration) 
{
    // TODO : where to put these parameters
    double Goal = 0.4;
    int NoMemoryTime = 10000;
    
    double Epsilon = DecreasingStepSize(Iteration, NoMemoryTime);
    double Denom;
    
    if(Ratio > Goal)
    {
        Denom = 1 - Goal;
    } 
    else
    {
        Denom = Goal; 
    }
    
    double NewVariance = GRV.GetVariance();
    NewVariance *= (1 + Epsilon * (Ratio - Goal) / Denom / 20);
    GRV.SetVariance( NewVariance );
}