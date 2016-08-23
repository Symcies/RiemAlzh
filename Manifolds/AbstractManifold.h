#ifndef _AbstractManifold_h
#define _AbstractManifold_h


#include <memory>
#include <string>
#include "../RandomVariables/GaussianRandomVariable.h"
#include "../RandomVariables/AbstractRandomVariable.h"
#include <map>


/// TODO : GetGeodesicDerivative, GetGeodesic and ComputeParallelTransport
/// ... shouldn't have std::vector<double> W0 as argument
/// ... W0 should be a realization !

class AbstractManifold {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    typedef std::map< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariableMap;
    typedef std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariable;
    typedef std::map<std::string, std::vector<double>> Realizations;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    AbstractManifold();
    ~AbstractManifold();


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    inline const double GetDimension() { return m_Dimension; }

    /*
    virtual inline const std::vector<double> GetGeodesicDerivative(double TimePoint, const std::shared_ptr<Realizations>& R) = 0;

    virtual inline const std::vector<double> GetGeodesic(double TimePoint, const std::shared_ptr<Realizations>& R) = 0;
    */

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compute the geodesic
    template <typename TimePoint>
    virtual const std::vector<double> ComputeGeodesic(std::vector<double> P0, double T0, std::vector<double> V0, TimePoint T) = 0;

    /// Compute the geodesic derivative
    template <typename TimePoint>
    virtual const std::vector<double> ComputeGeodesicDerivative(std::vector<double> P0, double T0, std::vector<double> V0, TimePoint T) = 0;

    /// Compute the parallel curve
    template <typename TimePoint>
    virtual const std::vector<double> ComputeParallelCurveBIS(std::vector<double> P0, double T0, std::vector<double> V0, std::vector<double> W, TimePoint T) = 0;


    /// Compute the parallel transport
    virtual std::vector<double> ComputeParallelTransport(double T, std::vector<double> W0, const std::shared_ptr<Realizations>& R ) = 0;

    /// Compute the parallel curve ;which is the Riemannian Exponential of the parallel Transport
    /// TO DO : To be removed
    virtual std::vector<double> ComputeParallelCurve(double TimePoint, std::vector<double> W0, const std::shared_ptr<Realizations>& R) = 0;

    /// Get any vector transformation  wrt the metric (used in the householder method)
    virtual std::vector<double> ComputeMetricTransformation(std::vector<double> VectorToTransform, std::vector<double> ApplicationPoint) = 0;

    /// Compute the scalar product corresponding to the manifold metric
    virtual double ComputeScalarProduct(std::vector<double> U, std::vector<double> V, std::vector<double> ApplicationPoint) = 0;

protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Dimension of the Riemanian Manifold
    unsigned int m_Dimension;

};


#endif //_AbstractManifold_h
