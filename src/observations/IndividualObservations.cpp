#include "IndividualObservations.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

IndividualObservations
::IndividualObservations(VectorType time_points)
{
  time_points_ = time_points;
  time_points_num_ = time_points_.size();

  landmarks_presence_ = false;
  cog_score_presence_ = false;
}


IndividualObservations::~IndividualObservations()
{

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void IndividualObservations::AddCognitiveScores(std::vector<VectorType> cog_scores)
{
  cog_scores_ = cog_scores;
  cog_score_presence_ = true;
}


void IndividualObservations::AddLandmarks(std::vector<VectorType> Landmarks) 
{
  landmarks_ = Landmarks;
  landmarks_presence_ = true;
}
