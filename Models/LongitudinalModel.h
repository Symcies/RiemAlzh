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

    LongitudinalModel(unsigned int NbIndependentComponents);
    ~LongitudinalModel();


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Change a realization in the model or manifold when it is not a pointer to the realization in the algorithm
    virtual void SetRealization(std::string Name, double Realization);

    /// Get the candidate random variable corresponding to a realization
    virtual std::shared_ptr< AbstractRandomVariable > GetCandidateRandomVariable(const std::string Name, const double Realization);

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the random variables : Population-wide and subject-specific
    virtual void InitializeRandomVariables();

    /// Initialize the model parameters, if any, and the manifold parameters, if any
    /// In the propagation case, only the Delta(k) (Manifold param) are initialized
    virtual void InitializeModelParameters(const std::shared_ptr< Realizations >& R);

     /// Update the sufficient statistics according to the model variables / parameters 
    virtual SufficientStatisticsVector GetSufficientStatistics(const std::shared_ptr<Realizations>& R, const std::shared_ptr<Data>& D);

    /// Update the fixed effects thanks to the approximation step of the algorithm
    virtual void UpdateRandomVariables(const SufficientStatisticsVector& SufficientStatistics, const std::shared_ptr<Data>& D);

    /// Compute the likelihood of the model
    virtual double ComputeLikelihood(const std::shared_ptr<Realizations>& R, const std::shared_ptr<Data>& D);

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

    /// Compute Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns>
    void ComputeOrthonormalBasis( const std::shared_ptr<Realizations>& R); // TODO : Use a library to do it faster

    /// Compute the A Matrix used to get the space shifts
    void ComputeAMatrix( const std::shared_ptr<Realizations>& R); // TODO : Use a library to do it faster

    // Compute the space shifts
    void ComputeSpaceShifts(const std::shared_ptr<Realizations>& R, const int NumberOfSubjects); // TODO : Use a library to do it faster


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

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


};


#endif //_LongitudinalModel_h
