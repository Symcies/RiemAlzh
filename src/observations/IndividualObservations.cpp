#include "IndividualObservations.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

IndividualObservations
::IndividualObservations(VectorType TimePoints) 
{
  m_TimePoints = TimePoints;
  m_NumberOfTimePoints = m_TimePoints.size();
  
  m_LandmarksPresence = false;
  m_CognitiveScoresPresence = false;
}


IndividualObservations
::~IndividualObservations() 
{
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
IndividualObservations
::AddCognitiveScores(std::vector<VectorType> CognitiveScores) 
{
  m_CognitiveScores = CognitiveScores;
  m_CognitiveScoresPresence = true;
}


void
IndividualObservations
::AddLandmarks(std::vector<VectorType> Landmarks) 
{
  m_Landmarks = Landmarks;
  m_LandmarksPresence = true;
}