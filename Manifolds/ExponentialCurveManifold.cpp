#include "ExponentialCurveManifold.h"

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
::ComputeParallelCurve(VectorType P0, double T0, VectorType V0, VectorType SpaceShift,
                       double TimePoint, VectorType Delta) 
{
    /// Initialization
    VectorType ParallelCurve(SpaceShift.size());
    double InitialPosition = P0(0);
    double InitialVelocity = V0(0);
    auto IterS = SpaceShift.begin(), IterD = Delta.begin(), IterP = ParallelCurve.begin();
    
    /// Compute coordinates
    for(    ; IterS != SpaceShift.end() && IterD != Delta.end() && IterP != ParallelCurve.end(); ++IterS, ++IterD, ++IterP)
    {
        double Val = *IterS / (InitialPosition * exp(*IterD)) + *IterD - TimePoint / InitialPosition;
        // TODO : BE FUCKING CAREFUL ABOUT " - TimePoint" OR " + TimePoint" !!
        *IterP = InitialPosition * exp(Val);
        //std::cout << *IterS / (InitialPosition * exp(*IterD)) << " & " << *IterD << " & " << - TimePoint / InitialPosition << std::endl;
        //std::cout << "W_i & 1/(P0*exp) : " << *IterS << " & " << 1 / (InitialPosition * exp(*IterD)) << std::endl;
        //std::cout << *IterP << std::endl;
    }
    
    return ParallelCurve;
}


ExponentialCurveManifold::VectorType
ExponentialCurveManifold
::ComputeParallelCurve(VectorType P0, double T0, VectorType V0, VectorType SpaceShift,
                       double TimePoint) 
{
    throw std::invalid_argument( "Wrong overloaded function ComputeParallelCurve called - in Exponential Curve Manifold" );
}

ExponentialCurveManifold::VectorType
ExponentialCurveManifold
::GetVelocityTransformToEuclideanSpace(VectorType P0, double T0, VectorType V0, VectorType Delta) 
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
::GetVelocityTransformToEuclideanSpace(VectorType P0, double T0, VectorType V0) 
{
    throw std::invalid_argument( "Wrong overloaded function GetVelocityTransformToEuclidianSpace called - in ExponentialCurve Manifold" );
}

double 
ExponentialCurveManifold
::ComputeScalarProduct(VectorType U, VectorType V, VectorType ApplicationPoint) 
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
