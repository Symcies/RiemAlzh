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
  
  unsigned int GetNumberOfTimePoints() const { return m_NumberOfTimePoints; }
  
  ScalarType GetTimePoint(unsigned int TimePointNumber) const { return m_TimePoints(TimePointNumber); }
  
  VectorType GetTimePoints() const { return m_TimePoints; }
  
  const VectorType& GetLandmark(unsigned int TimePointNumber) const { return m_Landmarks.at(TimePointNumber); }
  
  const VectorType& GetCognitiveScore(unsigned int TimePointNumber) const { return m_CognitiveScores.at(TimePointNumber); }
  
  const bool LandmarksPresence() const { return m_LandmarksPresence; }
  
  const bool CognitiveScoresPresence() const { return m_CognitiveScoresPresence; }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Add the cognitive scores
  void AddCognitiveScores(std::vector<VectorType> CognitiveScores);

  /// Add the landmarks
  void AddLandmarks(std::vector<VectorType> Landmarks) ;
  
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
  
  /// Number of timePoints
  unsigned int m_NumberOfTimePoints;
  
  /// List of patient observations
  VectorType m_TimePoints;
  
  /// List of cognitive scores - listed according to the observations
  std::vector<VectorType> m_CognitiveScores;
  
  /// Presence of cognitive scores
  bool m_CognitiveScoresPresence;

  /// List of landmarks - listed according to the observations
  std::vector<VectorType> m_Landmarks;
  
  /// Presence of landmarks
  bool m_LandmarksPresence;
  
};

