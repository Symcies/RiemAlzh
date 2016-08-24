
#ifndef _PropagationManifold_h
#define _PropagationManifold_h


#include "AbstractManifold.h"


class PropagationManifold : public AbstractManifold {
public:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    PropagationManifold(unsigned int NumberDimension, std::shared_ptr<AbstractBaseManifold>& BM);


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compute the geodesic
    virtual std::vector<double> ComputeGeodesic(std::vector<double> P0, double T0, std::vector<double> V0, double TimePoint);

    /// Compute the geodesic derivative
    virtual std::vector<double> ComputeGeodesicDerivative(std::vector<double> P0, double T0, std::vector<double> V0, double TimePoint);

    /// Compute the parallel curve
    virtual std::vector<double> ComputeParallelCurve(std::vector<double> P0, double T0, std::vector<double> V0,
                                                     std::vector<double> SpaceShift, double TimePoint );

    /// Compute the parallel transport
    virtual std::vector<double> ComputeParallelTransport(std::vector<double> P0, double T0, std::vector<double> V0,
                                                         std::vector<double> SpaceShift, double TimePoint);

    /// Get V0 transformation  wrt the metric at the application point P0 (used in the householder method)
    virtual std::vector<double> GetVelocityTransformToEuclideanSpace(std::vector<double> P0, double T0, std::vector<double> V0);

    ////////// TODO : CHECK WHERE TO PUT THE FOLLOWING FUNCTIONS
    /// Compute the scalar product corresponding to the manifold metric
    virtual double ComputeScalarProduct(std::vector<double> U, std::vector<double> V, std::vector<double> ApplicationPoint);


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<double> GetPropagationParameters();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

};


#endif //_PropagationManifold_h
