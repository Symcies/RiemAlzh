#pragma once

#include "AbstractBaseManifold.h"

class ExponentialBaseManifold : public AbstractBaseManifold {
public:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    ExponentialBaseManifold();
    virtual ~ExponentialBaseManifold();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
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
    ExponentialBaseManifold& operator=(const ExponentialBaseManifold&);
    ExponentialBaseManifold(const ExponentialBaseManifold&);

};
