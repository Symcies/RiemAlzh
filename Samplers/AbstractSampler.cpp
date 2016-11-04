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
    
    return 1.0 / pow(Epsilon, 0.6);
}

void
AbstractSampler
::UpdatePropositionDistributionVariance(GaussianRandomVariable &GRV, double Ratio, int Iteration) 
{
    // TODO : where to put these parameters
    double Goal = 0.3;
    double MovePossibility = 30.0; // Maximum percentage of move far from CurrentVariance 
    int NoMemoryTime = 400;
    
    double CurrentVariance = GRV.GetVariance();
    double Descent = DecreasingStepSize(Iteration, NoMemoryTime);
    double NewVariance = 1.0 + MovePossibility * (Goal - Ratio) * Descent * CurrentVariance;
    NewVariance *= CurrentVariance;
    GRV.SetVariance( NewVariance );
}