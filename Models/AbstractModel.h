#ifndef _AbstractModel_h
#define _AbstractModel_h

#include <memory>

#include "../RandomVariables/AbstractRandomVariable.h"
#include "../RandomVariables/LaplaceRandomVariable.h"
#include "../Manifolds/AbstractManifold.h"
#include <algorithm>

class AbstractModel {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    typedef std::vector< std::vector< std::pair< std::vector<double>, double> > > Data;
    typedef std::vector< std::pair< std::vector<double>, double> > IndividualData;
    typedef std::map< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariableMap;
    typedef std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariable;
    typedef std::map<std::string, std::vector<double>> Realizations;
    typedef std::vector< std::vector< double >> SufficientStatisticsVector;


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    inline void SetManifold(std::shared_ptr< AbstractManifold>& M) { m_Manifold = M; }

    RandomVariable GetRandomVariable(std::string name);

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the random variables : Population-wide and subject-specific
    virtual void InitializeRandomVariables() = 0;

    /// Initialize the model parameters, if any, and the manifold parameters, if any
    virtual void InitializeModelParameters(const std::shared_ptr< Realizations >& R) = 0;

    /// Update the sufficient statistics according to the model variables / parameters 
    virtual SufficientStatisticsVector GetSufficientStatistics(const std::shared_ptr<Realizations>& R, const std::shared_ptr<Data>& D) = 0;

    /// Update the fixed effects thanks to the approximation step of the algorithm
    virtual void UpdateRandomVariables(const std::vector< std::vector< double >>& SufficientStatistics, const std::shared_ptr<Data>& D) = 0;

    /// Compute the likelihood of the model
    virtual double ComputeLikelihood(const std::shared_ptr<Realizations>& R, const std::shared_ptr<Data>& D) = 0;

    /// Simulate data according to the model
    virtual Data* SimulateData(int NumberOfSubjects, int MinObs, int MaxObs) = 0;

    /// Simulate some random variable realizations
    Realizations* SimulateRealizations(int NumberOfSubjects);


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the true parameters to simulate data according to it - these parameters are unknown to the algo
    virtual void InitializeFakeRandomVariables() = 0;


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Riemanian manifold
    std::shared_ptr< AbstractManifold > m_Manifold;

    /// Random variables shared among the population
    RandomVariableMap m_PopulationRandomVariables;

    /// Random variables that are subject-specific
    RandomVariableMap m_IndividualRandomVariables;

    /// Random variables that are related to the manifold
    RandomVariableMap m_ManifoldRandomVariables;


};


#endif //_AbstractModel_h
