#include "IndividualObservations.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

IndividualObservations::IndividualObservations(VectorType time_points, int id)
{
  time_points_ = time_points;
  time_points_num_ = time_points_.size();

  id_ = id;

  landmarks_presence_ = false;
  cog_score_presence_ = false;
}


IndividualObservations::IndividualObservations(const IndividualObservations& obs)
{
  id_                 = obs.id_;
  rid_                = obs.rid_;
  time_points_num_    = obs.time_points_num_;
  time_points_        = obs.time_points_;
  cog_scores_         = obs.cog_scores_;
  cog_score_presence_ = obs.cog_score_presence_;
  landmarks_          = obs.landmarks_;
  landmarks_presence_ = obs.landmarks_presence_;
}

IndividualObservations::~IndividualObservations()
{

}

IndividualObservations& IndividualObservations::operator=(const IndividualObservations & obs){
    id_                 = obs.id_;
    rid_                = obs.rid_;
    time_points_num_    = obs.time_points_num_;
    time_points_        = obs.time_points_;
    cog_scores_         = obs.cog_scores_;
    cog_score_presence_ = obs.cog_score_presence_;
    landmarks_          = obs.landmarks_;
    landmarks_presence_ = obs.landmarks_presence_;
    return *this ;
};
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
