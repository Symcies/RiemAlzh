#ifndef _TestModel_h
#define _TestModel_h

#include "AbstractModel.h"

class TestModel : public AbstractModel {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    TestModel();
    ~TestModel();
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
  
    /// Initialize the model
    virtual void Initialize(const std::shared_ptr<const Data> D);
    
    /// Update the parameters
    virtual void UpdateParameters(const std::shared_ptr<Realizations>& R, const std::vector<std::string> Names = {"All"});
    

    /// Update the sufficient statistics according to the model variables / parameters 
    virtual SufficientStatisticsVector GetSufficientStatistics(const std::shared_ptr<Realizations>& R, const std::shared_ptr<Data>& D);

    /// Update the fixed effects thanks to the approximation step of the algorithm
    virtual void UpdateRandomVariables(const SufficientStatisticsVector& StochSufficientStatistics, const std::shared_ptr<Data>& D);
    
    
    /// Compute the log likelihood of the model
    virtual double ComputeLogLikelihood(const std::shared_ptr<Realizations>& R, const std::shared_ptr<Data>& D);
    
    /// Compute the log likelihood of the model for a particular individual
    virtual double ComputeIndividualLogLikelihood(const std::shared_ptr<Realizations>& R, 
                                                  const std::shared_ptr<Data>& D, const int SubjectNumber);
    
    /// Simulate data according to the model
    virtual Data SimulateData(int NumberOfSubjects, int MinObs, int MaxObs);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Outputs
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Compute Outputs
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

    /// Get a coefficient of subject i
    double GetA(int i, const std::shared_ptr<Realizations>& R);
    
    /// Get b coefficient of subject i
    double GetB(int i, const std::shared_ptr<Realizations>& R);
    
    /// get c coefficient (population ranom variable)
    double GetC(const std::shared_ptr<Realizations> R);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Noise associated to the model
    std::shared_ptr<GaussianRandomVariable> m_Noise;
    
};


#endif //_TestModel_h
