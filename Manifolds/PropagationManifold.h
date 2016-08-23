
#ifndef _PropagationManifold_h
#define _PropagationManifold_h


#include "AbstractManifold.h"


class PropagationManifold : public AbstractManifold {
public:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    typedef std::vector<std::shared_ptr< AbstractRandomVariable >> PointersVector;


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    PropagationManifold(unsigned int NumberDimension);


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    virtual inline const std::vector<double> GetGeodesicDerivative(double TimePoint, const std::shared_ptr<Realizations>& R);

    virtual inline const std::vector<double> GetGeodesic(double TimePoint, const std::shared_ptr<Realizations>& R);



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compute the parallel transport
    virtual std::vector<double> ComputeParallelTransport(double T, std::vector<double> W0, const std::shared_ptr<Realizations>& R );

    /// Compute the parallel transport
    virtual std::vector<double> ComputeParallelCurve(double TimePoint, std::vector<double> W0, const std::shared_ptr<Realizations>& R);

    /// Get any vector transformation  wrt the metric (used in the householder method)
    virtual std::vector<double> ComputeMetricTransformation(std::vector<double> VectorToTransform, std::vector<double> ApplicationPoint);

    /// Compute the scalar product corresponding to the manifold metric
    virtual double ComputeScalarProduct(std::vector<double> U, std::vector<double> V, std::vector<double> ApplicationPoint);


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compute the one dimensional geodesic
    double ComputeOneDimensionalGeodesic(double P0, double T0, double V0, double T);

    /// Compute the one dimensional geodesic derivative
    double ComputeOneDimensionalGeodesicDerivative(double P0, double T0, double V0, double T);

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

};


#endif //_PropagationManifold_h
