#ifndef _AbstractManifold_h
#define _AbstractManifold_h


#include <memory>
#include <string>
#include "../RandomVariables/GaussianRandomVariable.h"
#include "../RandomVariables/AbstractRandomVariable.h"
#include <map>

class AbstractManifold {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    typedef std::map< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariableMap;
    typedef std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariable;


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    AbstractManifold();
    ~AbstractManifold();


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    inline const double GetDimension() { return m_Dimension; }

    inline RandomVariableMap GetManifoldRandomVariables() { return m_ManifoldRandomVariables; }

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the random variables related to the manifold
    virtual void InitializeRandomVariables() = 0;

    /// Compute the parallel curve
    virtual std::vector<double> ComputeParallelCurve(double TimePoint, std::vector<double> W0, std::map<std::string, double> Realization) = 0;

protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Dimension of the Riemanian Manifold
    unsigned int m_Dimension;


    /// Random variables related to the manifold to be sent to the model
    RandomVariableMap m_ManifoldRandomVariables;
};


#endif //_AbstractManifold_h
