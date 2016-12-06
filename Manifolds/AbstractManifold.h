#ifndef _AbstractManifold_h
#define _AbstractManifold_h


typedef double ScalarType;

#include <memory>
#include <string>
#include "../Tests/TestAssert.h"
#include "../RandomVariables/GaussianRandomVariable.h"
#include "../RandomVariables/AbstractRandomVariable.h"
#include "BaseManifold/AbstractBaseManifold.h"
#include "../LinearAlgebra/LinearAlgebra.h"
#include <map>




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
    virtual VectorType ComputeGeodesic(VectorType P0, double T0, VectorType V0, double TimePoint);

    /// Compute the geodesic derivative
    virtual VectorType ComputeGeodesicDerivative(VectorType P0, double T0, VectorType V0, double TimePoint);

    /// Compute the parallel curve
    virtual VectorType ComputeParallelCurve(VectorType P0, double T0, VectorType V0,
                                                     VectorType SpaceShift, double TimePoint );

    /// Get V0 transformation  wrt the metric at the application point P0 (used in the householder method)
    virtual VectorType GetVelocityTransformToEuclideanSpace(VectorType P0, double T0, VectorType V0) = 0;

    /// Compute the scalar product corresponding to the manifold metric
    virtual double ComputeScalarProduct(VectorType U, VectorType V, VectorType ApplicationPoint) = 0;

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
