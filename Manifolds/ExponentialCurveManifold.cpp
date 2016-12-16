#include "ExponentialCurveManifold.h"
#include <cassert>
////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////



ExponentialCurveManifold
::ExponentialCurveManifold(unsigned int NumberOfDimension) 
{
    m_Dimension = NumberOfDimension;   
}


ExponentialCurveManifold
::~ExponentialCurveManifold() 
{
    
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

ExponentialCurveManifold::VectorType
ExponentialCurveManifold
::ComputeParallelCurve(VectorType& P0, double T0, VectorType& V0, VectorType& SpaceShift,
                       double TimePoint, VectorType& Delta) 
{
    
    assert(SpaceShift.size() == Delta.size());
    
    /// Initialization
    VectorType ParallelCurve(SpaceShift.size());
    double InitialPosition = P0(0);
    ScalarType * p = ParallelCurve.memptr();
    ScalarType * s = SpaceShift.memptr();
    ScalarType * d = Delta.memptr();
    auto N = ParallelCurve.size();
    
#pragma omp simd
    for(size_t i = 0; i < N; ++i)
        p[i] = InitialPosition * exp(s[i] / (InitialPosition * exp(d[i])) + d[i] - TimePoint / InitialPosition);
    
    return ParallelCurve;
}


ExponentialCurveManifold::VectorType
ExponentialCurveManifold
::ComputeParallelCurve(VectorType& P0, double T0, VectorType& V0, VectorType& SpaceShift,
                       double TimePoint) 
{
    throw std::invalid_argument( "Wrong overloaded function ComputeParallelCurve called - in Exponential Curve Manifold" );
}

ExponentialCurveManifold::VectorType
ExponentialCurveManifold
::GetVelocityTransformToEuclideanSpace(VectorType& P0, double T0, VectorType& V0, VectorType& Delta) 
{
    
    /// Initialization
    VectorType TransformedVelocity(Delta.size());
    double InitialPosition = P0(0);
    double InitialVelocity = V0(0);
    double Ratio = InitialVelocity  / (InitialPosition * InitialPosition);
    
    /// Compute the coordinates
    auto IterTV = TransformedVelocity.begin(), IterD = Delta.begin();
    for(    ; IterTV != TransformedVelocity.end() && IterD != Delta.end(); ++IterTV, ++IterD)
    {
        *IterTV = Ratio * exp( - *IterD );
    }
    
    return TransformedVelocity;
}


ExponentialCurveManifold::VectorType
ExponentialCurveManifold
::GetVelocityTransformToEuclideanSpace(VectorType& P0, double T0, VectorType& V0) 
{
    throw std::invalid_argument( "Wrong overloaded function GetVelocityTransformToEuclidianSpace called - in ExponentialCurve Manifold" );
}

double 
ExponentialCurveManifold
::ComputeScalarProduct(VectorType& U, VectorType& V, VectorType& ApplicationPoint) 
{
    TestAssert::WarningEquality_Object(U.size(), V.size(), "ExponentialCurveManifold > ComputeScalarProduct");
    TestAssert::WarningEquality_Object(U.size(), ApplicationPoint.size(), "ExponentialCurveManifold > ComputeScalarProduct");
    
    double ScalarProduct = 0.0;
    
    auto IterU = U.begin(), IterV = V.begin(), IterA = ApplicationPoint.begin();
    for(    ; IterU != U.end() && IterV != V.end() && IterA !=ApplicationPoint.end()
            ; ++IterU, ++IterV, ++IterA)
    {
        ScalarProduct += *IterU * *IterV / (*IterA * *IterA);
    }
    
    return ScalarProduct;
}
