#include "AbstractSampler.h"

AbstractSampler::~AbstractSampler(){}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Methods :
////////////////////////////////////////////////////////////////////////////////////////////////////

double AbstractSampler::DecreasingStepSize(int iter, int no_memory_time)
{
    double epsilon = std::max(1, iter - no_memory_time);
    double lambda = 0.51; // lambda should belong to ]1/2 ; 1]
    return 1.0 / pow(epsilon, lambda);
}
