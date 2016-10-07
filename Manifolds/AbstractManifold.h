#ifndef _AbstractManifold_h
#define _AbstractManifold_h


#include <memory>
#include <string>
#include "../RandomVariables/GaussianRandomVariable.h"
#include "../RandomVariables/AbstractRandomVariable.h"
#include "BaseManifold/AbstractBaseManifold.h"
#include <map>




class AbstractManifold {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    typedef std::map<std::string, double> Parameters;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    AbstractManifold();
    ~AbstractManifold();


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    inline const double GetDimension() { return m_Dimension; }


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compute the geodesic
    virtual std::vector<double> ComputeGeodesic(std::vector<double> P0, double T0, std::vector<double> V0, double TimePoint);

    /// Compute the geodesic derivative
    virtual std::vector<double> ComputeGeodesicDerivative(std::vector<double> P0, double T0, std::vector<double> V0, double TimePoint);

    /// Compute the parallel transport
    virtual std::vector<double> ComputeParallelTransport(std::vector<double> P0, double T0, std::vector<double> V0,
                                                         std::vector<double> SpaceShift, double TimePoint);
    /// Compute the parallel curve
    virtual std::vector<double> ComputeParallelCurve(std::vector<double> P0, double T0, std::vector<double> V0,
                                                     std::vector<double> SpaceShift, double TimePoint );

    /// Get V0 transformation  wrt the metric at the application point P0 (used in the householder method)
    virtual std::vector<double> GetVelocityTransformToEuclideanSpace(std::vector<double> P0, double T0, std::vector<double> V0) = 0;

    /// Compute the scalar product corresponding to the manifold metric
    virtual double ComputeScalarProduct(std::vector<double> U, std::vector<double> V, std::vector<double> ApplicationPoint) = 0;

protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Dimension of the Riemanian Manifold
    unsigned int m_Dimension;

    /// Base Manifold
    std::shared_ptr<AbstractBaseManifold> m_BaseManifold;
    
};


#endif //_AbstractManifold_h
