#ifndef _BlockedGibbsSampler2_h
#define _BlockedGibbsSampler2_h

#include <vector>
#include <string>
#include <map>
#include "AbstractSampler.h"
#include "../Models/AbstractModel.h"
#include "../Parameters/CandidateRandomVariables.h"
#include "../RandomVariables/AbstractRandomVariable.h"


class BlockedGibbsSampler2 : public AbstractSampler {
public:
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    BlockedGibbsSampler2();


    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

   
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the sampler
    virtual void InitializeSampler(const std::shared_ptr<MultiRealizations>& R);
    
    // Sample a new variable thanks to the sampler
    // The model cannot be constant because we modify some of its parameters (m_Orthonormal Basis for instance)
    virtual MultiRealizations Sample(const std::shared_ptr<MultiRealizations>& R, std::shared_ptr<AbstractModel>& M,
                                     const std::shared_ptr<Data>& D, int IterationNumber);


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    /// Sample the population variables
    MultiRealizations SamplePopulation(const std::shared_ptr<MultiRealizations>& R, std::shared_ptr<AbstractModel>& M, const std::shared_ptr<Data>& D);
    
    /// Get current population realizations
    UniqueRealizations GetCurrentPopulationRealizations(const std::shared_ptr<MultiRealizations>& R);
    
    /// Get current realizations of individual i
    UniqueRealizations GetCurrentIndividualRealizations(const std::shared_ptr<MultiRealizations>& R, int i);
    
    /// Get Candidatepopulation realizations
    UniqueRealizations GetCandidateRealizations(const std::shared_ptr<UniqueRealizations>& R, int Number);
    
    
    /// Sample the individual variables
    MultiRealizations SampleIndividual(int i, const std::shared_ptr<MultiRealizations>& R, std::shared_ptr<AbstractModel>& M,
                         const std::shared_ptr<Data>& D);
    

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Names of the population variables
    std::vector<std::string> m_NamePopulationVariables;
    
    /// Names of the individual variables
    std::vector<std::string> m_NameIndividualVariables;
};


#endif //_BlockedGibbsSampler_h