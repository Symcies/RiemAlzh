#ifndef _MeshworkModel_h
#define _MeshworkModel_h

#include "AbstractModel.h"

class MeshworkModel : public AbstractModel{
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    MeshworkModel(io::ModelSettings& MS);
    ~MeshworkModel();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the model
    virtual void Initialize(const Observations& Obs);
    
    /// Initialize the variance of the proposition distribution
    virtual ScalarType InitializePropositionDistributionVariance(std::string Name) const;
        
    /// Update parameters ; some model-specifid private members need to be initilize, m_Orthogonal Basis for instance
    /// This update can depend on the parameter that has changed, provided by the Name argument
    virtual void UpdateModel(const Realizations& R, int Type, const std::vector<std::string> Names = {"All"});

    /// Update the sufficient statistics according to the model variables / parameters 
    virtual SufficientStatisticsVector GetSufficientStatistics(const Realizations& R, const Observations& Obs);
    
    /// Update the fixed effects thanks to the approximation step of the algorithm
    virtual void UpdateRandomVariables(const SufficientStatisticsVector& StochSufficientStatistics);
    
    /// Compute the log likelihood of the model
    /// Using the log likelihood may have computational reason - for instance when the likelihood is too small
    virtual ScalarType ComputeLogLikelihood(const Observations &Obs);
    
    /// Compute the log likelihood of the model for a particular individual
    virtual ScalarType ComputeIndividualLogLikelihood(const IndividualObservations& Obs, const int SubjectNumber);
    
    /// Simulate data according to the model
    virtual Observations SimulateData(io::DataSettings& DS);
    
    /// Define the sampler block used in the gibbs sampler (should it be here?)
    virtual std::vector<SamplerBlock> GetSamplerBlocks() const;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Outputs
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Compute Outputs
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
    void ComputeThicknesses(const Realizations& R);
        
    /// Compute Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns>
    void ComputeOrthonormalBasis(); 

    /// Compute the A Matrix used to get the space shifts
    void ComputeAMatrix(const Realizations& R);

    /// Compute the space shifts
    void ComputeSpaceShifts(const Realizations& R); 
    
    /// Compute the block p0 * exp(delta_k)
    void ComputeBlock();
    
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
    unsigned int m_NbIndependentSources;
    
    /// Kernel Matrix K
    MatrixType m_InvertKernelMatrix;
    
    /// Interpolation Matrix to calculate any interpolation
    MatrixType m_InterpolationMatrix;
    
    /// Initial position of the model P0 = exp(R.at("P0")(0))
    VectorType m_Thicknesses;
    
    /// Interpolation coefficients of delta
    VectorType m_Deltas;
    
    /// Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns> (Basis orthogonal to gamma0_deriv(T0)
    MatrixType orthog_basis_;
    
    /// A Matrix vec<A1, ..., A(N)> where Ai is vec<Ns> (Ai is a column)
    MatrixType a_matrix_;

    /// Space shifts w(i) of the model
    MatrixType space_shifts_;
    
    /// Real time of observation of each individual
    std::vector<VectorType> indiv_obs_date_;
    
    /// Time reparametrization of each individual
    std::vector<VectorType> indiv_time_points_;
        
    /// Block1 corresponds to p0 * exp(Delta)
    VectorType block1_;
    
   
    

};


#endif //_MeshworkModel_h
