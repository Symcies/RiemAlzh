#pragma once

#include <cmath>

class AbstractBaseManifold {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    AbstractBaseManifold(){};
    virtual ~AbstractBaseManifold(){};

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compute the geodesic
    virtual double ComputeGeodesic(double p0, double t0, double v0, double time_point) = 0;

    /// Compute the geodesic derivative
    virtual double ComputeGeodesicDerivative(double p0, double t0, double v0, double time_point) = 0;

    /// Compute the parallel curve
    virtual double ComputeParallelCurve(double p0, double t0, double v0, double space_shift, double time_point) = 0;

    /// Compute the scalar product
    virtual double ComputeScalarProduct(double u, double v, double application_point) = 0;

private:
    /// Copy constructor, private to prevent copy
    AbstractBaseManifold(const AbstractBaseManifold &);

    /// Assignment operator, private to prevent copy
    AbstractBaseManifold& operator=(const AbstractBaseManifold &);

};
