#ifndef _NetworkPropagationModel_h
#define _NetworkPropagationModel_h


#include "AbstractModel.h"

class NetworkPropagationModel : public AbstractModel {
public: 
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    typedef std::vector< std::vector<double>> ControlPoints;
    
    // TODO : TO BE CHANGED ABSOLUTELY
    typedef std::vector<std::vector< double >> Matrix;
    typedef double Data;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the model : The random variables, the Kernel matrix and the interpolation matrix
    virtual void Initialize(const std::shared_ptr<ControlPoints>& P, const std::shared_ptr<Data>& D);
    
    /// Update parameters of the model, if any has to be updated
    virtual void UpdateParameters(const std::shared_ptr<MultiRealizations>& R, std::string Name = "All");
    
    // TODO : TO BE CHANGED ABSOLUTELLY : this is not how the random variables are updated GENERALLY
    // TODO : In fact, this was made because it is not generic by now as for the algorithm maximization step
    /// Update the fixed effects thanks to the approximation step of the algorithm
    virtual void UpdateRandomVariables(const SufficientStatisticsVector& StochSufficientStatistics, const std::shared_ptr<Data>& D) = 0;

    /// Compute the likelihood of the model
    virtual double ComputeLikelihood(const std::shared_ptr<MultiRealizations>& R, const std::shared_ptr<Data>& D, 
                                     const std::pair<std::string, int> NameRandomVariable = std::pair<std::string, int> ("All", 0)) = 0;


    /// Compute the log likelihood of the model
    /// Using the log likelihood may have computational reason - for instance when the likelihood is too small
    virtual double ComputeLogLikelihood(const std::shared_ptr<MultiRealizations>& R, const std::shared_ptr<Data>& D, 
                                        const std::pair<std::string, int> NameRandomVariable = std::pair<std::string, int> ("All", 0)) = 0;
    
    /// Compute the log likelihood of the model for a particular individual
    virtual double ComputeIndividualLogLikelihood(const std::shared_ptr<MultiRealizations>& R, 
                                                  const std::shared_ptr<Data>& D, const int SubjectNumber) = 0;

    /// Simulate data according to the model
    virtual Data SimulateData(int NumberOfSubjects, int MinObs, int MaxObs) = 0;

    
    
protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Compute the outputs
    virtual void ComputeOutputs() = 0;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Interpolation Matrix to calculate any point translation
    Matrix m_InterpolationMatrix;
    
    /// Kernel Matrix K 
    Matrix m_KernelMatrix;
    
    /// Noise associated to the model
    std::shared_ptr<GaussianRandomVariable> m_Noise;
    
    
    
};


#endif //_NetworkPropagationModel_h
