#ifndef _LinearManifold_h
#define _LinearManifold_h

#include "AbstractManifold.h"

class LinearManifold : public AbstractManifold{
public:
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    LinearManifold(unsigned int NumberOfDimension);
    ~LinearManifold();


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Compute the parallel curve
    virtual VectorType ComputeParallelCurve(VectorType& P0, double T0, VectorType& V0, VectorType& SpaceShift, double TimePoint, VectorType& Delta);
    
    /// Compute the parallel curve  - NULL propagation
    virtual VectorType ComputeParallelCurve(VectorType& P0, double T0, VectorType& V0, VectorType& SpaceShift, double TimePoint);
    
    /// Get V0 transformation  wrt the metric at the application point P0 (used in the householder method)
    virtual VectorType GetVelocityTransformToEuclideanSpace(VectorType& P0, double T0, VectorType& V0, VectorType& Delta);
    
    /// Get V0 transformation  wrt the metric at the application point P0 (used in the householder method)  - NULL propagation
    virtual VectorType GetVelocityTransformToEuclideanSpace(VectorType& P0, double T0, VectorType& V0);

    /// Compute the scalar product corresponding to the manifold metric
    virtual double ComputeScalarProduct(VectorType& U, VectorType& V, VectorType& ApplicationPoint);
};


#endif //_LinearManifold_h
