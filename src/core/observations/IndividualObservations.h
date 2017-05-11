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

  IndividualObservations(VectorType time_points, int id);
  ~IndividualObservations();

  IndividualObservations(const IndividualObservations &);
  IndividualObservations& operator=(const IndividualObservations &);

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  unsigned int GetNumberOfTimePoints() const { return time_points_num_; }

  unsigned int GetId() const { return id_; }

  ScalarType GetTimePoint(unsigned int time_points_num) const { return time_points_(time_points_num); }

  VectorType GetTimePoints() const { return time_points_; }

  const VectorType& GetLandmark(unsigned int time_points_num) const { return landmarks_.at(time_points_num); }

  const VectorType& GetCognitiveScore(unsigned int time_points_num) const { return cog_scores_.at(time_points_num); }

  const bool LandmarksPresence() const { return landmarks_presence_; }

  const bool CognitiveScoresPresence() const { return cog_score_presence_; }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Add the cognitive scores
  void AddCognitiveScores(std::vector<VectorType> cog_scores);

  /// Add the landmarks
  void AddLandmarks(std::vector<VectorType> landmarks) ;

private:

  /// ID of the patient
  unsigned int id_;

  /// RID (ADNI feature) of the patient
  unsigned int rid_;

  /// Number of timePoints
  unsigned int time_points_num_;

  /// List of patient observations
  VectorType time_points_;

  /// List of cognitive scores - listed according to the observations
  std::vector<VectorType> cog_scores_;

  /// Presence of cognitive scores
  bool cog_score_presence_;

  /// List of landmarks - listed according to the observations
  std::vector<VectorType> landmarks_;

  /// Presence of landmarks
  bool landmarks_presence_;



};
