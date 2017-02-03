#ifndef _MeshworkModel_h
#define _MeshworkModel_h

#include "ModelSettings.h"
#include "AbstractModel.h"


class MeshworkModel : public AbstractModel {
public:
     ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    typedef typename LinearAlgebra<ScalarType>::MatrixType MatrixType;
    typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    MeshworkModel(const ModelSettings& MS);
    ~MeshworkModel();
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
  
        /// Initialize the model : random variables, interpolation matrix, parameters
    virtual void Initialize(const Data& D);
    
    /// Update the model parameters != random variables parameters
    virtual void UpdateModel(const Realizations &R, int Type, const std::vector<std::string> Names = {"All"});
    
    /// Simulate data according to the model and the parameters
    virtual Data SimulateData(int NumberOfSubjects, int MinObs, int  MaxObs);
    
    /// Compute the log likelihood of the model
    virtual double ComputeLogLikelihood(const Realizations& R, const Data& D);
    
    /// Compute the log likelihood of the model for a given subject
    virtual double ComputeIndividualLogLikelihood(const Realizations& R, const Data& D, const int SubjectNumber);

    /// Get the sufficient statistics of the model
    virtual SufficientStatisticsVector GetSufficientStatistics(const Realizations& R, const Data& D);
    
    /// Update the random variables <=> the parameters of the model
    virtual void UpdateRandomVariables(const SufficientStatisticsVector& SS, const Data& D);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Outputs
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Compute outputs
    virtual void DisplayOutputs(const Realizations& R);
    
    /// Save the data into a file
    virtual void SaveData(unsigned int IterationNumber, const Realizations& R);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
  
    /// Initialize the true parameters to simulate data according to it - these parameters are unknown to the algo
    virtual void InitializeFakeRandomVariables();
    
    
protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Noise model
    std::shared_ptr<GaussianRandomVariable> m_Noise;
    
    /// Dimension of the manifold
    unsigned int m_ManifoldDimension;
    
    /// Number of indenpendent sources
    unsigned int m_NbIndependentSources;
    
    /// Number of control points
    unsigned int m_NbControlPoints;
    
    /// Thickess values P(k) at each patch
    VectorType m_Thicknesses;
    
    /// Time translation values Delta(k) at each patch
    VectorType m_Deltas;

    /// Subject time point
    MatrixType m_SubjectTimePoints;
    
    /// Subjects observation date
    MatrixType m_SubjectObervationDate;
    
    /// Block 1 corresponding to p_k + exp(d_k)
    VectorType m_Block1;
    

    
    /// Kernel Matrix K
    MatrixType m_InvertKernelMatrix;
    
    /// Interpolation Matrix to calculate any interpolation
    MatrixType m_InterpolationMatrix;
    
    
    
    /// Space shifts of the model
    std::vector<VectorType> m_SpaceShifts;
    
    
};


#endif //_MeshworkModel_h
