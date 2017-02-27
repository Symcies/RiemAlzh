#ifndef _UnivariateModel_h
#define _UnivariateModel_h


#include "AbstractModel.h"
#include <tuple>

class UnivariateModel : public AbstractModel {
public:
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    UnivariateModel(std::shared_ptr<AbstractBaseManifold>& M);
    ~UnivariateModel();
    

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Initialize the model
    virtual void Initialize(const Data& D);
    
    
    /// Update parameters ; some model-specifid private members need to be initilize, m_Orthogonal Basis for instance
    /// This update can depend on the parameter that has changed, provided by the Name argument
    virtual void UpdateModel(const Real &R, int Type, const std::vector<std::string> Names = {"All"});

    /// Update the sufficient statistics according to the model variables / parameters 
    virtual SufficientStatisticsVector GetSufficientStatistics(const Real& R, 
                                                               const Data& D);

    // TODO : TO BE CHANGED ABSOLUTELLY : this is not how the random variables are updated GENERALLY
    // TODO : In fact, this was made because it is not generic by now as for the algorithm maximization step
    /// Update the fixed effects thanks to the approximation step of the algorithm
    virtual void UpdateRandomVariables(const SufficientStatisticsVector& StochSufficientStatistics, 
                                       const Data& D);
    
    /// Compute the log likelihood of the model
    /// Using the log likelihood may have computational reason, e.g. when the likelihood is too small
    virtual double ComputeLogLikelihood(const Data &D);

    /// Compute the log likelihood of the model for a particular individual
    double ComputeIndividualLogLikelihood(const Data &D, const int SubjectNumber);

    
    /// Simulate data according to the model
    virtual Data SimulateData(int NumberOfSubjects, int MinObs, int MaxObs);

    /// Save the data into a file
    virtual void SaveData(unsigned int IterationNumber, const Real& R);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Outputs
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Compute Outputs
    virtual void DisplayOutputs(const Real& R);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the true parameters to simulate data according to it - these parameters are unknown to the algo
    virtual void InitializeFakeRandomVariables();

private : 
    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Get the subject time point psi_i(t) = exp(ksi_i) * (t - T0 - tau_i) - T0
    std::function<double(double)> GetSubjectTimePoint(const int SubjectNumber, const Real& R);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Last calculated Likelihood - and the corresponding realizations
    /// Bool : if last calculation was generic. Double : last likelihood value. Realizations : last realizations
    std::tuple<bool, double, Real> m_LastLogLikelihood;
    
    /// Noise associated to the model
    std::shared_ptr<GaussianRandomVariable> m_Noise;
    
    /// Abstract Base Manifold - because do not need an AbstractManifold
    std::shared_ptr<AbstractBaseManifold> m_BaseManifold;
};


#endif //_UnivariateModel_h
