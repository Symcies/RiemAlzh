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
    double Lambda = 0.51; // Lambda should belong to ]1/2 ; 1]
    return 1.0 / pow(Epsilon, Lambda);
}

void
AbstractSampler
::UpdatePropositionDistributionVariance(GaussianRandomVariable &GRV, double Ratio, int Iteration) 
{
    // TODO : where to put these parameters
    double Goal = 0.259;
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
    
    double CurrentVariance = GRV.GetVariance();
    double NewVariance = CurrentVariance * (1 + Epsilon * (Ratio - Goal) / Denom / 5);
    NewVariance = std::max(NewVariance, 0.0000001);
    GRV.SetVariance( NewVariance );
}