#ifndef _LongitudinalModel_h
#define _LongitudinalModel_h


#include "AbstractModel.h"
#include "../Utilities/MatrixFunctions.h"

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

    /// Initialize the random variables : Population-wide and subject-specific
    virtual void InitializeRandomVariables();

     /// Update the sufficient statistics according to the model variables / parameters 
    virtual SufficientStatisticsVector GetSufficientStatistics(const std::shared_ptr<Realizations>& R, const std::shared_ptr<Data>& D);

    /// Update the fixed effects thanks to the approximation step of the algorithm
    virtual void UpdateRandomVariables(const SufficientStatisticsVector& SufficientStatistics, const std::shared_ptr<Data>& D);


    /// Compute the likelihood of the model
    virtual double ComputeLikelihood(const std::shared_ptr<Realizations>& R, const std::shared_ptr<Data>& D, const std::pair<std::string, int> NameRandomVariable);

    /// Compute the log likelihood of the model
    /// Using the log likelihood may have computational reason - for instance when the likelihood is too small
    virtual double ComputeLogLikelihood(const std::shared_ptr<Realizations>& R, const std::shared_ptr<Data>& D, const std::pair<std::string, int> NameRandomVariable);

    /// Simulate data according to the model
    virtual Data SimulateData(int NumberOfSubjects, int MinObs, int MaxObs);


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the true parameters to simulate data according to it - these parameters are unknown to the algo
    virtual void InitializeFakeRandomVariables();



protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize Population random variables
    void InitializePopulationRandomVariables();

    /// Initialize Individual random variables
    void InitializeIndividualRandomVariables();

    /// Initialize Manifold random variables
    void InitializeManifoldRandomVariables();

    // TODO : Peut -être faire des get pour toutes les différentes réalisations?
    /// Get the initial position = gamma(t0)
    std::vector<double> GetInitialPosition(const std::shared_ptr<Realizations>& R);

    /// Get the initial velocity = diff(gamma(t0))
    std::vector<double> GetInitialVelocity(const std::shared_ptr<Realizations>& R);

    /// Get the propagation = (delta(k))
    std::vector<double> GetPropagationCoefficients(const std::shared_ptr<Realizations>& R);

    /// Compute Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns>
    void ComputeOrthonormalBasis( const std::shared_ptr<Realizations>& R); // TODO : Use a library to do it faster

    /// Compute the A Matrix used to get the space shifts
    void ComputeAMatrix( const std::shared_ptr<Realizations>& R); // TODO : Use a library to do it faster

    // Compute the space shifts
    void ComputeSpaceShifts(const std::shared_ptr<Realizations>& R); // TODO : Use a library to do it faster

    /// Compute the Likelihood the most generic way, without simplification
    double ComputeLogLikelihoodGeneric(const std::shared_ptr<Realizations>& R, const std::shared_ptr<Data>& D);

    /// Compute the likelihood keeping the term of the specific individual
    double ComputeLogLikelihoodIndividual(const std::shared_ptr<Realizations>& R, const std::shared_ptr<Data>& D, const int SubjectNumber);

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Last calculated Likelihood - and the corresponding realizations
    /// Bool : if last calculation was generic. Double : last likelihood value. Realizations : last realizations
    std::tuple<bool, double, Realizations> m_LastLogLikelihood;

    /// Number of independent components
    unsigned int m_NbIndependentComponents;

    /// Noise model
    std::shared_ptr< GaussianRandomVariable > m_Noise;

    /// Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns> (Basis orthogonal to gamma0_deriv(T0)
    std::vector< std::vector< double >> m_OrthogonalBasis; // TODO : Initialize somewhere

    /// A Matrix vec<A1, ..., A(N)> where Ai is vec<Ns> (Ai is a column)
    std::vector< std::vector< double >> m_AMatrix; // TODO : Initialize somewhere

    /// Space shifts w(i) of the model
    std::map< std::string, std::vector< double >> m_SpaceShifts; // TODO : Initialize somewhere

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Output(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compute outputs
    virtual void ComputeOutputs();

};


#endif //_LongitudinalModel_h
