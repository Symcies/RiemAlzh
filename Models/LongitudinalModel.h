#ifndef _LongitudinalModel_h
#define _LongitudinalModel_h


#include "AbstractModel.h"
#include <functional>
#include <vector>

class LongitudinalModel : public AbstractModel {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    LongitudinalModel(const unsigned int NbIndependentComponents, std::shared_ptr<AbstractManifold>& M);
    ~LongitudinalModel();


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the model
    virtual void Initialize(const Data& D);
    
    
    /// Initialize parameters ; some model-specifid private members need to be initilize, m_Orthogonal Basis for instance
    virtual void UpdateModel(const Realizations &R, const std::vector<std::string> Names = {"All"});
    
     /// Update the sufficient statistics according to the model variables / parameters 
    virtual SufficientStatisticsVector GetSufficientStatistics(const Realizations& R, const Data& D);

    /// Update the fixed effects thanks to the approximation step of the algorithm
    virtual void UpdateRandomVariables(const SufficientStatisticsVector& SS, 
                                       const Data& D);
    
    /// Compute the log likelihood of the model
    /// Using the log likelihood may have computational reason - for instance when the likelihood is too small
    virtual double ComputeLogLikelihood(const Realizations& R, const Data& D);
    
    /// Compute the log likelihood of the model for a particular individual
    double ComputeIndividualLogLikelihood(const Realizations& R, const Data& D, const int SubjectNumber);
    
    /// Simulate data according to the model
    virtual Data SimulateData(int NumberOfSubjects, int MinObs, int MaxObs);


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the true parameters to simulate data according to it - these parameters are unknown to the algo
    virtual void InitializeFakeRandomVariables();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Outputs
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Compute Outputs
    virtual void DisplayOutputs();
    
    /// Save the data into a file
    virtual void SaveData(unsigned int IterationNumber, const Realizations &R);


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Get the initial time 
    double GetInitialTime();
    
    /// Get the propagation coefficients = (delta(k))
    VectorType GetPropagationCoefficients(const Realizations& R);
    
    /// Get the subject time point psi_i(t) = exp(ksi_i) * (t - T0 - tau_i) - T0
    std::function<double(double)> GetSubjectTimePoint(const int SubjectNumber, const Realizations& R);

    /// Compute Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns>
    void ComputeOrthonormalBasis( const Realizations& R);

    /// Compute the A Matrix used to get the space shifts
    void ComputeAMatrix( const Realizations& R);

    /// Compute the space shifts
    void ComputeSpaceShifts(const Realizations& R);
    
    // TODO : To delete and add to a manifold function, if possible
    /// Compute the geodesic
    VectorType ComputeGeodesic(double P0, double TimePoint, VectorType Delta, VectorType SpaceShift);
    
    /// Compute the geodesic transformation to get a euclidean scalar product
    VectorType ComputeGeodesicTransformation(double P0, double T0, double V0, VectorType Delta);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Manifold dimension
    unsigned int m_ManifoldDimension;
    
    /// Number of independent components
    unsigned int m_NbIndependentComponents;
    
    /// Sum of the observations - corresponds to the first sufficient statistic
    double m_SumObservations;
    
    /// Total number of observations throughout all the individuals
    double m_NbTotalOfObservations;
    
    /// Noise model
    std::shared_ptr< GaussianRandomVariable > m_Noise;

    /// Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns> (Basis orthogonal to gamma0_deriv(T0)
    std::vector< VectorType > m_OrthogonalBasis;

    /// A Matrix vec<A1, ..., A(N)> where Ai is vec<Ns> (Ai is a column)
    MatrixType m_AMatrix;

    /// Space shifts w(i) of the model
    std::map< std::string, VectorType> m_SpaceShifts;
    
};


#endif //_LongitudinalModel_h
