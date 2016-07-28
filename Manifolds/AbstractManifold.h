#ifndef _AbstractManifold_h
#define _AbstractManifold_h


#include <memory>
#include <string>
#include "../RandomVariables/GaussianRandomVariable.h"
#include "../RandomVariables/AbstractRandomVariable.h"
#include <map>


/// TODO : GetGeodesicDerivative, GetGeodesic and ComputeParallelTransport
/// ... shouldn't have std::vector<double> W0 as argument
/// ... W0 should be a realization !

class AbstractManifold {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    typedef std::map< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariableMap;
    typedef std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariable;
    typedef std::map<std::string, double> Realizations;

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

    virtual inline const std::vector<double> GetGeodesicDerivative(double TimePoint, const Realizations& R) = 0;

    virtual inline const std::vector<double> GetGeodesic(double TimePoint, const Realizations& R) = 0;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the random variables related to the manifold
    virtual void InitializeRandomVariables() = 0;

    /// Compute the parallel curve
    virtual std::vector<double> ComputeParallelCurve(double TimePoint, std::vector<double> W0, const Realizations& R) = 0;

    /// Get any vector transformation  wrt the metric (used in the householder method)
    virtual std::vector<double> ComputeMetricTransformation(std::vector<double> VectorToTransform, std::vector<double> ApplicationPoint) = 0;

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
