#ifndef _AbstractSampler_h
#define _AbstractSampler_h


#include <memory>
#include <string>
#include "../RandomVariables/GaussianRandomVariable.h"
#include "../RandomVariables/AbstractRandomVariable.h"
#include <map>
#include "../Models/AbstractModel.h"

class AbstractSampler {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    typedef std::vector< std::vector< std::pair< std::vector<double>, double> > > Data;
    typedef std::map< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariableMap;
    typedef std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariable;
    typedef std::map<std::string, std::vector<double>> Realizations;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    AbstractSampler();
    ~AbstractSampler();


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

   
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Sample a new variable thanks to the sampler
    virtual void Sample(RandomVariable& CurrentRV, std::shared_ptr< AbstractRandomVariable >& CandidateRV, std::shared_ptr<AbstractModel>& M, Realizations& R, const std::shared_ptr<Data>& D) = 0;


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////


};


#endif //_AbstractSampler_h
