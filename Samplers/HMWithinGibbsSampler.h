#ifndef _HMWithinGibbsSampler_h
#define _HMWithinGibbsSampler_h


#include <memory>
#include <string>
#include "../RandomVariables/GaussianRandomVariable.h"
#include "../RandomVariables/AbstractRandomVariable.h"
#include "AbstractSampler.h"
#include "../Models/AbstractModel.h"
#include <map>

class HMWithinGibbsSampler : public AbstractSampler{
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    typedef std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariable;
    typedef std::map<std::string, double> Realizations;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    HMWithinGibbsSampler();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Sample a new variable thanks to the sampler
    virtual void Sample(RandomVariable& CurrentRV, double& CurrentRealization, RandomVariable& CandidateRV, AbstractModel& M, Realizations& R, Data& D);
   
   
protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////


};


#endif //_HMWithinGibbsSampler_h
