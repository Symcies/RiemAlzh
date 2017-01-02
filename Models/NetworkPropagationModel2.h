#ifndef _NetworkPropagationModel2_h
#define _NetworkPropagationModel2_h


#include "AbstractModel.h"

class NetworkPropagationModel2 : public AbstractModel {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    typedef typename LinearAlgebra<ScalarType>::MatrixType MatrixType;
    typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    NetworkPropagationModel2(const unsigned int NbIndependentComponents, std::shared_ptr<AbstractManifold> M,
                             std::shared_ptr<MatrixType> KernelMatrix, std::shared_ptr<MatrixType> InterpolationMatrix);
    
    ~NetworkPropagationModel2();
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
  
    /// Initialize the model : random variables, interpolation matrix, parameters
    virtual void Initialize(const Data& D);
    
    /// Update the model parameters != random variables parameters
    virtual void UpdateModel(const Realizations &R,
                             const std::vector<std::string> Names = {"All"});
    
    /// Simulate data according to the model and the parameters
    virtual Data SimulateData(int NumberOfSubjects, int MinObs, int  MaxObs);
    
    /// Compute the log likelihood of the model
    virtual double ComputeLogLikelihood(const Realizations& R, const Data& D);
    
    /// Compute the log likelihood of the model for a given subject
    virtual double ComputeIndividualLogLikelihood(const Realizations& R, 
                                          const Data& D, const int SubjectNumber);

    /// Get the sufficient statistics of the model
    virtual SufficientStatisticsVector GetSufficientStatistics(const Realizations& R, const Data& D);
    
    /// Update the random variables <=> the parameters of the model
    virtual void UpdateRandomVariables(const SufficientStatisticsVector& SS, const Data& D);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Outputs
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Compute outputs
    virtual void DisplayOutputs();
    
    /// Save the data into a file
    virtual void SaveData(unsigned int IterationNumber);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
  
    /// Initialize the true parameters to simulate data according to it - these parameters are unknown to the algo
    virtual void InitializeFakeRandomVariables();
    
protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Get the propagation coefficients delta_tilde
    VectorType GetPropagationCoefficients(const Realizations& R);
    
    /// Get the timepoint reparametrization for a given subject
    std::function<double(double)> GetSubjectTimePoint(const int SubjectNumber, const Realizations& R);
    
    /// Compute the Interpolation coefficients
    void ComputeInterpolationCoefficients(const Realizations& R);
    
    /// Compute Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns>
    void ComputeOrthonormalBasis(const Realizations& R); 

    /// Compute the A Matrix used to get the space shifts
    void ComputeAMatrix(const Realizations& R); 

    /// Compute the space shifts
    void ComputeSpaceShifts(const Realizations& R);
 

    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Noise model
    std::shared_ptr< GaussianRandomVariable > m_Noise;
    
    /// Number of control points
    unsigned int m_NbControlPoints;
    
    /// Number of independent components
    unsigned int m_NbIndependentComponents;
    
    /// Sum of the observations - corresponds to the first sufficient statistic
    double m_SumObservations;
    
    /// Total number of observations throughout all the individuals
    double m_NbTotalOfObservations;
    
    /// Kernel Matrix K
    MatrixType m_InvertKernelMatrix;
    
    /// Interpolation Matrix to calculate any interpolation
    MatrixType m_InterpolationMatrix;
    
    /// Interpolation coefficients a_i
    VectorType m_InterpolationCoefficients;
    
    /// Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns> (Basis orthogonal to gamma0_deriv(T0)
    std::vector< VectorType > m_OrthogonalBasis;

    /// A Matrix vec<A1, ..., A(N)> where Ai is vec<Ns> (Ai is a column)
    MatrixType m_AMatrix;

    /// Space shifts w(i) of the model
    std::map< std::string, VectorType> m_SpaceShifts;
    
};


#endif //_NetworkPropagationModel2_h
