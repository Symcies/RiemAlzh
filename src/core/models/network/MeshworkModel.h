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

    MeshworkModel(io::ModelSettings& model_settings);
    ~MeshworkModel();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the model
    virtual void Initialize(const Observations& obs);

    /// Initialize the variance of the proposition distribution
    virtual ScalarType InitializePropositionDistributionVariance(std::string name) const;

    /// Update parameters ; some model-specifid private members need to be initilize, m_Orthogonal Basis for instance
    /// This update can depend on the parameter that has changed, provided by the name argument
    virtual void UpdateModel(const Realizations& reals, int type, const std::vector<std::string> names = {"All"});

    /// Update the sufficient statistics according to the model variables / parameters
    virtual SufficientStatisticsVector GetSufficientStatistics(const Realizations& reals, const Observations& obs);

    /// Update the fixed effects thanks to the approximation step of the algorithm
    virtual void UpdateRandomVariables(const SufficientStatisticsVector& stoch_sufficient_stats);

    /// Compute the log likelihood of the model
    /// Using the log likelihood may have computational reason - for instance when the likelihood is too small
    virtual ScalarType ComputeLogLikelihood(const Observations &obs);

    /// Compute the log likelihood of the model for a particular individual
    virtual ScalarType ComputeIndividualLogLikelihood(const IndividualObservations& obs, const int indiv_num);

    /// Simulate data according to the model
    virtual Observations SimulateData(io::DataSettings& data_settings);

    /// Define the sampler block used in the gibbs sampler (should it be here?)
    virtual std::vector<SamplerBlock> GetSamplerBlocks() const;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Outputs
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compute Outputs
    virtual void DisplayOutputs(const Realizations& reals);

    /// Save the data into a file
    virtual void SaveData(unsigned int IterationNumber, const Realizations& reals);

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
    void ComputeSubjectTimePoint(const Realizations& reals, const int indiv_num = -1);

    /// Compute the interpolation coefficients delta
    void ComputeDeltas(const Realizations& reals);

    /// Compute the interpolation coefficients beta
    void ComputeThicknesses(const Realizations& reals);

    /// Compute Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns>
    void ComputeOrthonormalBasis();

    /// Compute the A Matrix used to get the space shifts
    void ComputeAMatrix(const Realizations& reals);

    /// Compute the space shifts
    void ComputeSpaceShifts(const Realizations& reals);

    /// Compute the block p0 * exp(delta_k)
    void ComputeBlock();

    /// Compute the parallel curve
    VectorType ComputeParallelCurve(int indiv_num, int obs_num);

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    /// Noise model
    std::shared_ptr< GaussianRandomVariable > noise_;

    /// Number of control points
    unsigned int control_points_nb_;

    /// Number of independent components
    unsigned int indep_components_nb_;

    /// Kernel Matrix K
    MatrixType invert_kernel_matrix_;

    /// Interpolation Matrix to calculate any interpolation
    MatrixType interpolation_matrix_;

    /// Initial position of the model P0 = exp(reals.at("P0")(0))
    VectorType thickenesses_;

    /// Interpolation coefficients of delta
    VectorType deltas_;

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


private:
    MeshworkModel(const MeshworkModel &);
    MeshworkModel& operator=(const MeshworkModel &);


};


#endif //_MeshworkModel_h
