#pragma once

#include "LinearAlgebra.h"

typedef double ScalarType;

class IndividualObservations {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    typedef typename LinearAlgebra<ScalarType>::MatrixType MatrixType;
    typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    IndividualObservations(VectorType TimePoints);
    ~IndividualObservations();
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Add the cognitive scores
    void AddCognitiveScores(std::vector<VectorType> CognitiveScores) { m_CognitiveScores = CognitiveScores; }

    /// Add the landmarks
    void AddLandmarks(std::vector<VectorType> Landmarks) { m_Landmarks = Landmarks; }
    
protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// ID of the patient
    unsigned int m_ID;
    
    /// RID (ADNI feature) of the patient
    unsigned int m_RID;
    
    /// List of patient observations
    VectorType m_TimePoints;
    
    /// List of cognitive scores - listed according to the observations
    std::vector<VectorType> m_CognitiveScores;
    
    /// List of landmarks - listed according to the observations
    std::vector<VectorType> m_Landmarks;
    
};

