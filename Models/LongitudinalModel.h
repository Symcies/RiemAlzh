#ifndef _LongitudinalModel_h
#define _LongitudinalModel_h


#include "AbstractModel.h"
#include "../Utilities/MatrixFunctions.h"
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
    virtual void Initialize();
    
    
    /// Initialize parameters ; some model-specifid private members need to be initilize, m_Orthogonal Basis for instance
    virtual void UpdateParameters(const std::shared_ptr<MultiRealizations>& R, const std::vector<std::string> Names = {"All"});

    /// Get the parameters of the model
    virtual std::map< std::string, double > GetParameters();
    
     /// Update the sufficient statistics according to the model variables / parameters 
    virtual SufficientStatisticsVector GetSufficientStatistics(const std::shared_ptr<MultiRealizations>& R, const std::shared_ptr<Data>& D);

    /// Update the fixed effects thanks to the approximation step of the algorithm
    virtual void UpdateRandomVariables(const SufficientStatisticsVector& StochSufficientStatistics, const std::shared_ptr<Data>& D);
    
    /// Compute the log likelihood of the model
    /// Using the log likelihood may have computational reason - for instance when the likelihood is too small
    virtual double ComputeLogLikelihood(const std::shared_ptr<MultiRealizations>& R, const std::shared_ptr<Data>& D);
    
    /// Compute the log likelihood of the model for a particular individual
    double ComputeIndividualLogLikelihood(const std::shared_ptr<MultiRealizations>& R, 
                                          const std::shared_ptr<Data>& D, const int SubjectNumber);
    
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
    virtual void ComputeOutputs();


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Get the initial time 
    double GetInitialTime();
    
    /// Get the initial position = gamma(t0)
    VectorType GetInitialPosition(const std::shared_ptr<MultiRealizations>& R);

    /// Get the initial velocity = diff(gamma(t0))
    VectorType GetInitialVelocity(const std::shared_ptr<MultiRealizations>& R);

    /// Get the propagation coefficients = (delta(k))
    VectorType GetPropagationCoefficients(const std::shared_ptr<MultiRealizations>& R);
    
    /// Get the subject time point psi_i(t) = exp(ksi_i) * (t - T0 - tau_i) - T0
    std::function<double(double)> GetSubjectTimePoint(const int SubjectNumber, const std::shared_ptr<MultiRealizations>& R);

    /// Compute Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns>
    void ComputeOrthonormalBasis( const std::shared_ptr<MultiRealizations>& R);

    /// Compute the A Matrix used to get the space shifts
    void ComputeAMatrix( const std::shared_ptr<MultiRealizations>& R);

    // Compute the space shifts
    void ComputeSpaceShifts(const std::shared_ptr<MultiRealizations>& R);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Last calculated Likelihood - and the corresponding realizations
    /// Bool : if last calculation was generic. Double : last likelihood value. Realizations : last realizations
    std::tuple<bool, double, MultiRealizations> m_LastLogLikelihood;

    /// Number of independent components
    unsigned int m_NbIndependentComponents;

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
