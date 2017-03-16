#ifndef _BlockedGibbsSampler_h
#define _BlockedGibbsSampler_h

#include "AbstractSampler.h"
#include <algorithm>

class BlockedGibbsSampler : public AbstractSampler {

public:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Constructor
    BlockedGibbsSampler();

    /// Constructor 2
    BlockedGibbsSampler(unsigned int memoryless_sampling_time, double expected_acceptance_ratio);

    /// Destructor
    ~BlockedGibbsSampler();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the sampler
    virtual void InitializeSampler(Realizations& reals, AbstractModel &models);

    /// Sample new realizations
    virtual void Sample(Realizations& reals, AbstractModel& models, const Observations& obs);


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Sample one block
    void OneBlockSample(int block_num, Realizations& reals, AbstractModel& models, const Observations& obs);

    /// Compute Prior part of the ratio while updating the realization
    ScalarType ComputePriorRatioAndUpdateRealizations(Realizations& reals, const AbstractModel& models, const MiniBlock& var);


    /// Compute likelihood based on the block type
    VectorType ComputeLogLikelihood(AbstractModel& models, const Observations& obs);

    /// Get previously computed log likelihood
    double GetPreviousLogLikelihood();

    /// Update the last log likelihood computed
    void UpdateLastLogLikelihood(VectorType& computed_log_likelihood);

        /// Update the random variable
    void UpdateBlockRandomVariable(double acceptance_ratio, const MiniBlock& var);

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Current Block Type
    int cur_block_type_;

    /// Current iteration Number
    int cur_iter_;



    /// Last likelihood computed : reused
    VectorType last_likelihood_computed_;

    /// Parameters updated in the block
    std::vector<std::string> cur_block_params_;
    std::vector<int> cur_block_params_bis_;

    /// Parameters to recover back in the realization
    std::vector<std::tuple<std::string, unsigned int, ScalarType >> recover_params_;

    /// Uniform distrubution
    std::uniform_real_distribution<double> uniform_distrib_;

};


#endif //_BlockedGibbsSampler_h
