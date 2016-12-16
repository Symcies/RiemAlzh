#include "LinearManifold.h"
#include <cassert>


LinearManifold
::LinearManifold(unsigned int NumberOfDimension) 
{
    m_Dimension = NumberOfDimension;
}

LinearManifold
::~LinearManifold() 
{
    
}



LinearManifold::VectorType
LinearManifold
::ComputeParallelCurve(VectorType &P0, double T0, VectorType &V0, VectorType &SpaceShift,
                       double TimePoint, VectorType& Delta) 
{
    assert(SpaceShift.size() == Delta.size());
    
    /// Initialization
    VectorType ParallelCurve(SpaceShift.size());
    double InitialPosition = P0(0);
    double InitialVelocity = V0(0);
    ScalarType * p = ParallelCurve.memptr();
    ScalarType * s = SpaceShift.memptr();
    ScalarType * d = Delta.memptr();
    auto N = ParallelCurve.size();
    
#pragma omp simd
    for(size_t i = 0; i < N; ++i)
        p[i] = InitialPosition + s[i] + d[i]*InitialPosition - TimePoint;
    
    return ParallelCurve;
}


LinearManifold::VectorType
LinearManifold
::ComputeParallelCurve(VectorType &P0, double T0, VectorType &V0, VectorType &SpaceShift,
                       double TimePoint) 
{
    throw std::invalid_argument( "Wrong overloaded function ComputeParallelCurve called - in Exponential Curve Manifold" );
}

LinearManifold::VectorType
LinearManifold
::GetVelocityTransformToEuclideanSpace(VectorType &P0, double T0, VectorType &V0,
                                       VectorType &Delta) 
{
    double InitialVelocity = V0(0);
    return VectorType(Delta.size(), InitialVelocity);
}

LinearManifold::VectorType
LinearManifold
::GetVelocityTransformToEuclideanSpace(VectorType &P0, double T0, VectorType &V0) 
{
    throw std::invalid_argument( "Wrong overloaded function GetVelocityTransformToEuclidianSpace called - in ExponentialCurve Manifold" );
}


double 
LinearManifold
::ComputeScalarProduct(VectorType &U, VectorType &V, VectorType &ApplicationPoint) 
{
    double ScalarProduct = 0.0;
    
    auto IterU = U.begin(), IterV = V.begin();
    for(    ; IterU != U.end() && IterV != V.end()
            ; ++IterU, ++IterV)
    {
        ScalarProduct += *IterU * *IterV;
    }
    
    return ScalarProduct;
}