#ifndef _FastNetworkModel_h
#define _FastNetworkModel_h


#include "AbstractModel.h"
#include "ExponentialCurveManifold.h"

class FastNetworkModel : public AbstractModel {
public:
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    typedef typename LinearAlgebra<ScalarType>::MatrixType MatrixType;
    typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    FastNetworkModel(io::ModelSettings& MS);
    FastNetworkModel(const unsigned int NbIndependentComponents, std::shared_ptr<MatrixType> KernelMatrix, std::shared_ptr<MatrixType> InterpolationMatrix);
    ~FastNetworkModel();
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
  
    /// Initialize the model : random variables, interpolation matrix, parameters
    virtual void Initialize(const Observations& Obs);
    
        /// Initialize the variance of the proposition distribution
    virtual ScalarType InitializePropositionDistributionVariance(std::string Name) const;
    
    /// Update the model parameters != random variables parameters
    virtual void UpdateModel(const Realizations& R, int Type, const std::vector<std::string> Names = {"All"});
    
    /// Simulate data according to the model and the parameters
    virtual OldData SimulateData(io::DataSettings& DS);
    
    /// Compute the log likelihood of the model
    virtual double ComputeLogLikelihood(const OldData& D);
    
    /// Compute the log likelihood of the model for a given subject
    virtual double ComputeIndividualLogLikelihood(const IndividualObservations& Obs, const int SubjectNumber);

    /// Get the sufficient statistics of the model
    virtual SufficientStatisticsVector GetSufficientStatistics(const Realizations& R, const Observations& Obs);
    
    /// Update the random variables <=> the parameters of the model
    virtual void UpdateRandomVariables(const SufficientStatisticsVector& SS);
    
    /// Define the sampler block used in the gibbs sampler (should it be here?)
    virtual std::vector<SamplerBlock> GetSamplerBlocks() const ;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Outputs
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Compute outputs
    virtual void DisplayOutputs(const Realizations& R);
    
    /// Save the data into a file
    virtual void SaveData(unsigned int IterationNumber, const Realizations& R);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
  
    /// Initialize the true parameters to simulate data according to it - these parameters are unknown to the algo
    virtual void InitializeFakeRandomVariables();

protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Compute the subjects time points
    void ComputeSubjectTimePoint(const Realizations& R, const int SubjectNumber = -1);
    
    /// Compute the interpolation coefficients delta
    void ComputeDeltas(const Realizations& R);
    
    /// Compute the interpolation coefficients beta
    void ComputeNus(const Realizations& R);
        
    /// Compute Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns>
    void ComputeOrthonormalBasis(); 

    /// Compute the A Matrix used to get the space shifts
    void ComputeAMatrix(const Realizations& R);

    /// Compute the space shifts
    void ComputeSpaceShifts(const Realizations& R); 
    
    /// Compute the time reparametrizations
    
    /// Compute the block p0 * exp(delta_k)
    void ComputeBlock1();
    
    /// Compute the block nu_k / p0
    void ComputeBlock2();
    
    /// Compute the parallel curve
    VectorType ComputeParallelCurve(int SubjectNumber, int ObservationNumber);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
        
    /// Noise model
    std::shared_ptr< GaussianRandomVariable > m_Noise;
    
    /// Number of control points
    unsigned int m_NbControlPoints;
    
    /// Number of independent components
    unsigned int m_NbIndependentComponents;
        
    /// Kernel Matrix K
    MatrixType m_InvertKernelMatrix;
    
    /// Interpolation Matrix to calculate any interpolation
    MatrixType m_InterpolationMatrix;
    
    /// Initial position of the model P0 = exp(R.at("P0")(0))
    double m_P0;
    
    /// Interpolation coefficients of delta
    VectorType m_Deltas;
    
    /// Interpolation coefficients of beta
    VectorType m_Nus;
    
    /// Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns> (Basis orthogonal to gamma0_deriv(T0)
    MatrixType m_OrthogonalBasis;
    
    /// A Matrix vec<A1, ..., A(N)> where Ai is vec<Ns> (Ai is a column)
    MatrixType m_AMatrix;

    /// Space shifts w(i) of the model
    MatrixType m_SpaceShifts;
    
    /// Real time of observation of each individual
    std::vector<VectorType> m_IndividualObservationDate;
    
    /// Time reparametrization of each individual
    std::vector<VectorType> m_SubjectTimePoints;
        
    /// Block1 corresponds to p0 * exp(Delta)
    VectorType m_Block1;
    
    /// Block2 corresponds to vu_k / p0
    VectorType m_Block2;
    
    

};


#endif //_FastNetworkModel_h
