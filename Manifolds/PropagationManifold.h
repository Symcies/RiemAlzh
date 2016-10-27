
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

    // Compute the geodesic
    virtual VectorType ComputeGeodesic(VectorType P0, double T0, VectorType V0, double TimePoint, VectorType Delta);

    /// Compute the geodesic - NULL propagation
    virtual VectorType ComputeGeodesic(VectorType P0, double T0, VectorType V0, double TimePoint);

    /// Compute the geodesic derivative
    virtual VectorType ComputeGeodesicDerivative(VectorType P0, double T0, VectorType V0,
                                                          double TimePoint, VectorType Delta);

    /// Compute the geodesic derivative - NULL propagation
    virtual VectorType ComputeGeodesicDerivative(VectorType P0, double T0, VectorType V0, double TimePoint);

    /// Compute the parallel curve
    virtual VectorType ComputeParallelCurve(VectorType P0, double T0, VectorType V0,
                                                     VectorType SpaceShift, double TimePoint, VectorType Delta  );

    /// Compute the parallel curve - NULL propagation
    virtual VectorType ComputeParallelCurve(VectorType P0, double T0, VectorType V0,
                                                     VectorType SpaceShift, double TimePoint);

    /// Compute the parallel transport
    virtual VectorType ComputeParallelTransport(VectorType P0, double T0, VectorType V0,
                                                         VectorType SpaceShift, double TimePoint, VectorType Delta);

    /// Compute the parallel transport - NULL propagation
    virtual VectorType ComputeParallelTransport(VectorType P0, double T0, VectorType V0,
                                                         VectorType SpaceShift, double TimePoint);

    /// Get V0 transformation  wrt the metric at the application point P0 (used in the householder method) - NULL propagation
    virtual VectorType GetVelocityTransformToEuclideanSpace(VectorType P0, double T0, VectorType V0, VectorType Delta);

    /// Get V0 transformation  wrt the metric at the application point P0 (used in the householder method)
    virtual VectorType GetVelocityTransformToEuclideanSpace(VectorType P0, double T0, VectorType V0);

    /// Compute the scalar product corresponding to the manifold metric
    virtual double ComputeScalarProduct(VectorType U, VectorType V, VectorType ApplicationPoint);


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

};


#endif //_PropagationManifold_h
