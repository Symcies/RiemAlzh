#ifndef _NetworkPropagationModel_h
#define _NetworkPropagationModel_h

#include "AbstractModel.h"


class NetworkPropagationModel : public AbstractModel {
public:
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    typedef typename LinearAlgebra<ScalarType>::MatrixType MatrixType;
    typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    NetworkPropagationModel(const unsigned int NbIndependentComponents, std::shared_ptr<MatrixType> KernelMatrix, std::shared_ptr<MatrixType> InterpolationMatrix);
    ~NetworkPropagationModel();
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
  
    /// Initialize the model : random variables, interpolation matrix, parameters
    virtual void Initialize(const std::shared_ptr<const Data> D);
    
    /// Update the model parameters != random variables parameters
    virtual void UpdateParameters(const Realizations& R, 
                                  const std::vector<std::string> Names = {"All"});
    
    /// Simulate data according to the model and the parameters
    virtual Data SimulateData(int NumberOfSubjects, int MinObs, int  MaxObs);
    
    /// Compute the log likelihood of the model
    virtual double ComputeLogLikelihood(const std::shared_ptr<Realizations> R, 
                                        const std::shared_ptr<const Data> D);
    
    /// Compute the log likelihood of the model for a given subject
    virtual double ComputeIndividualLogLikelihood(const std::shared_ptr<Realizations> R, 
                                                  const std::shared_ptr<const Data> D, 
                                                  const int SubjectNumber);

    /// Get the sufficient statistics of the model
    virtual SufficientStatisticsVector GetSufficientStatistics(const std::shared_ptr<Realizations> R, 
                                                               const std::shared_ptr<const Data> D);
    
    /// Update the random variables <=> the parameters of the model
    virtual void UpdateRandomVariables(const SufficientStatisticsVector& SS, 
                                       const std::shared_ptr<const Data> D);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Outputs
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Compute outputs
    virtual void ComputeOutputs();
    
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
    
    /// Get the delta for all nodes
    VectorType GetDelta(const std::shared_ptr<Realizations> R);
    
    /// Get the beta - v0 related - for all nodes
    VectorType GetNu(const std::shared_ptr<Realizations> R);
    
     /// Get the timepoint reparametrization for a given subject
    std::function<double(double)> GetSubjectTimePoint(const int SubjectNumber, 
                                                      const std::shared_ptr<Realizations> R);
    
    /// Compute the interpolation coefficients delta
    void ComputeInterpoCoeffDelta(const std::shared_ptr<Realizations> R);
    
    /// Compute the interpolation coefficients beta
    void ComputeInterpoCoeffNu(const std::shared_ptr<Realizations> R);
        
    /// Compute Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns>
    void ComputeOrthonormalBasis( const std::shared_ptr<Realizations> R); 

    /// Compute the A Matrix used to get the space shifts
    void ComputeAMatrix( const std::shared_ptr<Realizations> R);

    /// Compute the space shifts
    void ComputeSpaceShifts(const std::shared_ptr<Realizations> R); 
    
    
    /// Compute the parallel curve
    VectorType ComputeParallelCurve(VectorType& P0, VectorType& Delta, VectorType& Nu, VectorType& SpaceShift, double Timepoint);
    
    
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
    
    /// Interpolation coefficients of delta
    VectorType m_InterpolationCoeffDelta;
    
    /// Interpolation coefficients of beta
    VectorType m_InterpolationCoeffNu;
    
    
    /// Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns> (Basis orthogonal to gamma0_deriv(T0)
    std::vector< VectorType > m_OrthogonalBasis;

    /// A Matrix vec<A1, ..., A(N)> where Ai is vec<Ns> (Ai is a column)
    MatrixType m_AMatrix;

    /// Space shifts w(i) of the model
    std::map< std::string, VectorType> m_SpaceShifts;
    
    
    double m_ManifoldDimension = 1827;
};


#endif //_NetworkPropagationModel_h
