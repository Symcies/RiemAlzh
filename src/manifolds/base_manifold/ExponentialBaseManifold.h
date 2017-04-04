#pragma once

#include "AbstractBaseManifold.h"

class ExponentialBaseManifold : public AbstractBaseManifold {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    ExponentialBaseManifold();
    virtual ~ExponentialBaseManifold();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compute the geodesic
    virtual double ComputeGeodesic(double p0, double t0, double v0, double time_point);

    /// Compute the geodesic derivative
    virtual double ComputeGeodesicDerivative(double p0, double t0, double v0, double time_point);

    /// Compute the parallel curve
    virtual double ComputeParallelCurve(double p0, double t0, double v0, double space_shift, double time_point);

    /// Compute the scalar product
    virtual double ComputeScalarProduct(double u, double v, double application_point);

private:
    /// Assignment operator, private to prevent copy
    ExponentialBaseManifold& operator=(const ExponentialBaseManifold&);

    /// Copy constructor, private to prevent copy
    ExponentialBaseManifold(const ExponentialBaseManifold&);

};
