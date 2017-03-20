#pragma once

#include "AbstractManifold.h"

class LinearManifold : public AbstractManifold{
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    LinearManifold(unsigned int dim_num);
    ~LinearManifold();


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compute the parallel curve
    virtual VectorType ComputeParallelCurve(VectorType& p0, double t0, VectorType& v0, VectorType& space_shift, double time_point, VectorType& delta);

    /// Compute the parallel curve  - NULL propagation
    virtual VectorType ComputeParallelCurve(VectorType& p0, double t0, VectorType& v0, VectorType& space_shift, double time_point);

    /// Get v0 transformation  wrt the metric at the application point p0 (used in the householder method)
    virtual VectorType GetVelocityTransformToEuclideanSpace(VectorType& p0, double t0, VectorType& v0, VectorType& delta);

    /// Get v0 transformation  wrt the metric at the application point p0 (used in the householder method)  - NULL propagation
    virtual VectorType GetVelocityTransformToEuclideanSpace(VectorType& p0, double t0, VectorType& v0);

    /// Compute the scalar product corresponding to the manifold metric
    virtual double ComputeScalarProduct(VectorType& u, VectorType& v, VectorType& application_point);
};
