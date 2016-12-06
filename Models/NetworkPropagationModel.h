#ifndef _NetworkPropagationModel_h
#define _NetworkPropagationModel_h

typedef double ScalarType;


#include "AbstractModel.h"
#include "../Manifolds/PropagationManifold.h"
#include "../LinearAlgebra/LinearAlgebra.h"

class NetworkPropagationModel : public AbstractModel {
public: 
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    typedef std::vector< VectorType > ControlPoints;
    typedef std::vector< VectorType > MeshPoints;
    typedef typename LinearAlgebra<ScalarType>::MatrixType MatrixType;
    typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
        
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    NetworkPropagationModel(const unsigned int NbIndependentComponents, std::shared_ptr<AbstractManifold>& M, 
                            std::shared_ptr<MatrixType>& KernelMatrix, std::shared_ptr<MatrixType>& InterpolationMatrix);
    ~NetworkPropagationModel();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the model : The random variables, the Kernel matrix and the interpolation matrix
    virtual void Initialize(const std::shared_ptr<Data> D);
    
    /// Update parameters of the model, if any has to be updated
    virtual void UpdateParameters(const std::shared_ptr<MultiRealizations>& R, const std::vector<std::string> Names = {"All"});
    
    /// Get the parameters of the model
    virtual std::map< std::string, double > GetParameters();
    
    /// Update the sufficient statistics according to the model variables / parameters 
    virtual SufficientStatisticsVector GetSufficientStatistics(const std::shared_ptr<MultiRealizations>& R, const std::shared_ptr<Data>& D);

    
    // TODO : TO BE CHANGED ABSOLUTELLY : this is not how the random variables are updated GENERALLY
    // TODO : In fact, this was made because it is not generic by now as for the algorithm maximization step
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
    // Output(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compute outputs
    virtual void ComputeOutputs();
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the true parameters to simulate data according to it - these parameters are unknown to the algo
    virtual void InitializeFakeRandomVariables();
    
    
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
    void ComputeOrthonormalBasis( const std::shared_ptr<MultiRealizations>& R); // TODO : Use a library to do it faster

    /// Compute the A Matrix used to get the space shifts
    void ComputeAMatrix( const std::shared_ptr<MultiRealizations>& R); // TODO : Use a library to do it faster

    /// Compute the space shifts
    void ComputeSpaceShifts(const std::shared_ptr<MultiRealizations>& R); // TODO : Use a library to do it faster
    
    /// Compute the interpolation coefficients
    void ComputeInterpolationCoefficients(const std::shared_ptr<MultiRealizations>& R);
    
 

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
    
    /// Interpolation Matrix to calculate any point translation
    MatrixType m_InterpolationMatrix;
    
    /// Interpolation coefficients a(i)
    VectorType m_InterpolationCoefficients;
    
    /// Kernel Matrix K 
    MatrixType m_InvertKernelMatrix;
    
    
};


#endif //_NetworkPropagationModel_h
