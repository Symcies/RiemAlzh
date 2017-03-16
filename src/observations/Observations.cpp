#include "Observations.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

Observations::Observations()
{
  tot_indiv_num_ = 0;
}


Observations::~Observations()
{

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void Observations::AddIndividualData(IndividualObservations& indiv_obs)
{
  data_.push_back(indiv_obs);
  indiv_obs_.push_back(indiv_obs.GetTimePoints());
}


void Observations::InitializeGlobalAttributes() 
{
  tot_indiv_num_ = data_.size();
  tot_obs_num_ = 0;
  cog_scores_tot_sum_ = 0;
  landmarks_tot_sum_ = 0;

  for(auto it = data_.begin(); it != data_.end(); ++it)
  {
    tot_obs_num_ += it->GetNumberOfTimePoints();
    for(size_t i = 0; i < it->GetNumberOfTimePoints(); ++i)
    {
      if(it->LandmarksPresence())
      {
        landmarks_tot_sum_ += it->GetLandmark(i).squared_magnitude();
      }
      if(it->CognitiveScoresPresence())
      {
        cog_scores_tot_sum_ += it->GetCognitiveScore(i).squared_magnitude();
      }
    }
  }
}
