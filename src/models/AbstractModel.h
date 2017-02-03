#ifndef _AbstractModel_h
#define _AbstractModel_h

#include <algorithm>
#include <iostream>
#include <fstream>
#include <functional>
#include <math.h>
#include <memory>
#include <unordered_map>

#include "ReadData.h"
#include "LinearAlgebra.h"
#include "MatrixFunctions.h"
#include "TestAssert.h"
#include "AbstractRandomVariable.h"
#include "LaplaceRandomVariable.h"
#include "AbstractManifold.h"


class AbstractModel {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

        
    typedef typename LinearAlgebra<ScalarType>::MatrixType MatrixType;
    typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;

    
    typedef std::vector< std::vector< std::pair< VectorType, double> > > Data;
    typedef std::vector< std::pair< VectorType, double> > IndividualData;
    //typedef std::map<std::string, VectorType> Realizations;
    typedef std::unordered_map<std::string, VectorType> Realizations;
    typedef std::vector<VectorType> SufficientStatisticsVector;


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    std::shared_ptr< AbstractRandomVariable > GetRandomVariable(std::string name);

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the model
    virtual void Initialize(const Data& D) = 0;
        
    /// Update parameters ; some model-specifid private members need to be initilize, m_Orthogonal Basis for instance
    /// This update can depend on the parameter that has changed, provided by the Name argument
    virtual void UpdateModel(const Realizations &R, int Type,
                             const std::vector<std::string> Names = {"All"}) = 0;

    /// Update the sufficient statistics according to the model variables / parameters 
    virtual SufficientStatisticsVector GetSufficientStatistics(const Realizations& R, const Data& D) = 0;
    
    /// Update the fixed effects thanks to the approximation step of the algorithm
    virtual void UpdateRandomVariables(const SufficientStatisticsVector& StochSufficientStatistics, 
                                       const Data& D) = 0;
    
    
    /// Compute the log likelihood of the model
    /// Using the log likelihood may have computational reason - for instance when the likelihood is too small
    virtual double ComputeLogLikelihood(const Realizations& R, const Data& D)= 0;
    
    /// Compute the log likelihood of the model for a particular individual
    virtual double ComputeIndividualLogLikelihood(const Realizations& R, const Data& D, 
                                                  const int SubjectNumber) = 0;
    
    /// Simulate data according to the model
    virtual Data SimulateData(int NumberOfSubjects, int MinObs, int MaxObs) = 0;

    /// Simulate some random variable realizations
    Realizations SimulateRealizations(int NumberOfSubjects);

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Outputs
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Compute Outputs
    virtual void DisplayOutputs(const Realizations& R) = 0;
    
    /// Save the data into a file
    virtual void SaveData(unsigned int IterationNumber, const Realizations& R) = 0;
    
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
    std::map< std::string, std::shared_ptr< AbstractRandomVariable >>  m_PopulationRandomVariables;

    /// Random variables that are subject-specific
    std::map< std::string, std::shared_ptr< AbstractRandomVariable >> m_IndividualRandomVariables;
    
    /// Output file
    std::ofstream m_OutputParameters;


};


#endif //_AbstractModel_h
