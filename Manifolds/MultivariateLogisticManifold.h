#ifndef _MultivariateLogisticManifold_h
#define _MultivariateLogisticManifold_h


#include <stdexcept>            // throw
#include <memory>
#include <GaussianRandomVariable.h>
#include "AbstractManifold.h"

class MultivariateLogisticManifold : public AbstractManifold {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    MultivariateLogisticManifold(int DimensionNumber, int NumberOfIndependentComponents);


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Getter(s) and Setter(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Get the number of independent components Ns
    inline int GetNumberOfIndependentComponents() const { return m_NumberOfIndependentComponents; }

    // Set the propagation coefficient
    void SetPropagationCoefficient(std::vector< std::shared_ptr< GaussianRandomVariable> > PropagationCoefficient);


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Compute the Riemannian metric
    double ComputeMetric(std::vector<double> u,std::vector<double> v, std::vector<double> p);

    // Compute parallel Curve
    std::vector<double> ComputeParallelCurve(double P0, double T0, double V0, std::vector<double> W0, double T);

    // Compute geodesic
    std::vector<double> ComputeGeodesic(double P0, double T0, double V0, double T);

    // Compute double Geodesic Derivative
    std::vector<double> ComputeGeodesicDerivative(double P0, double T0, double V0, double T);


protected:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Number of Ns statistically independent components
    int m_NumberOfIndependentComponents;

    // Pointer to the propagation coefficient
    std::vector< std::shared_ptr< GaussianRandomVariable >> m_PropagationCoefficient;


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Compute one dimensional geodesic
    double ComputeOneDimensionalGeodesic(double P0, double T0, double V0, double T);

    // Compute parallel transport of vector W0
    std::vector<double> ComputeParallelTransport(double P0, double T0, double V0, std::vector<double> W0, double T) ;

    // Compute one dimensionsional geodesic derivative
    double ComputeOneDimensionalGeodesicDerivative(double P0, double T0, double V0, double T);


};


#endif //_MultivariateLogisticManifold_h
