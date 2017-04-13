#pragma once

typedef double ScalarType;

#include <memory>

#include "AbstractBaseManifold.h"
#include "LinearAlgebra.h"


class AbstractManifold {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    typedef typename LinearAlgebra<ScalarType>::MatrixType MatrixType;
    typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    AbstractManifold(){};
    virtual ~AbstractManifold(){};

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    inline const double GetDimension() { return dimension_; }

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compute the geodesic
    virtual VectorType ComputeGeodesic(VectorType& p0, double t0,
                                       VectorType& v0, double time_point);

    /// Compute the geodesic derivative
    virtual VectorType ComputeGeodesicDerivative(VectorType& p0, double t0,
                                                VectorType& v0, double time_point);

    /// Compute the parallel curve
    virtual VectorType ComputeParallelCurve(VectorType& p0, double t0, VectorType& v0,
                                            VectorType& space_shift, double time_point );

    /// Get v0 transformation  wrt the metric at the application point p0 (used in the householder method)
    virtual VectorType GetVelocityTransformToEuclideanSpace(VectorType& p0,
                                                            double t0,
                                                            VectorType& v0) = 0;

    /// Compute the scalar product corresponding to the manifold metric
    virtual double ComputeScalarProduct(VectorType& u, VectorType& v,
                                        VectorType& application_point) = 0;

protected:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Dimension of the Riemanian Manifold
    unsigned int dimension_;

    /// Base Manifold
    std::shared_ptr<AbstractBaseManifold> base_manifold_;

private:
    /// Copy constructor, private to prevent copy
    AbstractManifold(const AbstractManifold &);

    /// Assignment operator, private to prevent copy
    AbstractManifold& operator=(const AbstractManifold &);
};