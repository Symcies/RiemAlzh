#include "Observations.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

Observations
::Observations() 
{
  m_NumberOfSubjects = 0;
}


Observations
::~Observations() 
{
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void
Observations
::AddIndividualData(IndividualObservations& ID) 
{
  m_Data.push_back(ID);
  m_IndividualObservations.push_back(ID.GetTimePoints());
}


void 
Observations
::InitializeGlobalAttributes() 
{
  m_NumberOfSubjects = m_Data.size();
  m_TotalNumberOfObservations = 0;
  m_TotalSumOfCognitiveScores = 0;
  m_TotalSumOfLandmarks = 0;
  
  for(auto it = m_Data.begin(); it != m_Data.end(); ++it)
  {
    m_TotalNumberOfObservations += it->GetNumberOfTimePoints();
    for(size_t i = 0; i < it->GetNumberOfTimePoints(); ++i)
    {
      if(it->LandmarksPresence())      
      {
        m_TotalSumOfLandmarks += it->GetLandmark(i).squared_magnitude();
      }
      if(it->CognitiveScoresPresence()) 
      {
        m_TotalSumOfCognitiveScores += it->GetCognitiveScore(i).squared_magnitude();
      }
    }
  }
}