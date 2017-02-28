#pragma once

typedef double ScalarType;

#include "LinearAlgebra.h"
#include "IndividualObservations.h"

class Observations {
public:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    typedef typename LinearAlgebra<ScalarType>::MatrixType MatrixType;
    typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    Observations();
    ~Observations();
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    unsigned int GetNumberOfSubjects() const { return m_NumberOfSubjects; }
    
    ScalarType GetTotalNumberOfObservations() const { return m_TotalNumberOfObservations; }
    
    ScalarType GetTotalSumOfCognitiveScores() const { return m_TotalSumOfCognitiveScores; }
    
    ScalarType GetTotalSumOfLandmarks() const { return m_TotalSumOfLandmarks; }
    
    std::vector<VectorType> GetObservations() const { return m_IndividualObservations; }
    
    unsigned int GetNumberOfTimePoints(unsigned int SubjectNumber) const { return m_Data.at(SubjectNumber).GetNumberOfTimePoints(); }
    
    ScalarType GetSubjectTimePoint(unsigned int SubjectNumber, unsigned int TimePointNumber) const { 
        return m_Data.at(SubjectNumber).GetTimePoint(TimePointNumber); 
    }
    
    const VectorType& GetSubjectLandmark(unsigned int SubjectNumber, unsigned int TimePointNumber) const {
        return m_Data.at(SubjectNumber).GetLandmark(TimePointNumber);
    }
    
    const VectorType& GetSubjectCognitiveScore(unsigned int SubjectNumber, unsigned int TimePointNumber) const {
        return m_Data.at(SubjectNumber).GetCognitiveScore(TimePointNumber);
    }
    
    const IndividualObservations& GetSubjectObservations(unsigned int SubjectNumber) const { return m_Data.at(SubjectNumber); }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Add the observations of a new patient
    void AddIndividualData(IndividualObservations& ID);
    
    /// Initialize the cross-subjects attributes
    void InitializeGlobalAttributes();

protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Subject attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Individuals
    std::vector<IndividualObservations> m_Data;
    
    /// Individual time points
    std::vector<VectorType> m_IndividualObservations;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Cross-subject attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Total number of subjects
    unsigned int m_NumberOfSubjects;
    
    /// Total number of observations
    ScalarType m_TotalNumberOfObservations;
    
    /// Sum of the cognitive scores
    ScalarType m_TotalSumOfCognitiveScores;
    
    /// Sum of the landmarks
    ScalarType m_TotalSumOfLandmarks;
    
    
};
