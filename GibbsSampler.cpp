#include "GibbsSampler.h"




double ProposalDistribution(double p_prop, double p_init, double var) {
    // Hasting-Metropolis Random Walk : in the case of a normal distribution
    return 1.0/sqrt(var * 2 * M_PI) * exp( -1.0/2.0 * (p_prop - p_init) * (p_prop - p_init));
}

double AprioriDistribution()