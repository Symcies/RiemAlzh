#pragma once

typedef double ScalarType;

#include <fstream>
#include <iostream>
#include <string>

#include "tinyxml2.h"

namespace io {

class DataSettings {
public:
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////

  DataSettings(const char *XMLFile);

  ~DataSettings();

  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////

  std::string GetPathToGroup() const { return group_path_; }

  std::string GetPathToTimepoints() const { return timepoints_path_; }

  std::string GetPathToCognitiveScores() const { return cog_scores_path_; }

  unsigned int GetCognitiveScoresDimension() const { return cog_scores_dim_; }

  std::string GetPathToLandmarks() const { return landmarks_path; }

  unsigned int GetLandmarksDimension() const { return landmarks_dim_; }

  unsigned int GetNumberOfSimulatedSubjects() const { return subjects_total_num_; }

  unsigned int GetMinimumNumberOfObservations() const { return min_observation_num_; }

  unsigned int GetMaximumNumberOfObservations() const { return max_observation_num_; }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////

  bool IsReal() const { return is_data_real_; }

  bool LandmarkPresence() const { return are_landmarks_present_; }

  bool CognitiveScoresPresence() const { return are_cog_scores_present_; }


private:

  /// Read real data if true
  bool is_data_real_;

  /// Presence of landmarks
  bool are_landmarks_present_;

  /// Presence of cognitive scores
  bool are_cog_scores_present_;

  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// Attribute(s) for real data :
  ////////////////////////////////////////////////////////////////////////////////////////////////

  /// Path to data
  std::string data_path_;

  /// Path to group csv
  std::string group_path_;

  /// Path to timepoint csv
  std::string timepoints_path_;

  /// Path to cognitive scores (csv file)
  std::string cog_scores_path_;

  /// Dimension of the cognitive scores
  unsigned int cog_scores_dim_;

  /// Path to landmarks (csv file)
  std::string landmarks_path;

  /// Dimension of the landmarks
  unsigned int landmarks_dim_;


  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// Attribute(s) for simulated data :
  ////////////////////////////////////////////////////////////////////////////////////////////////

  /// NumberOfIndividuals to simulate
  unsigned int subjects_total_num_;

  /// Minimum number of observations
  unsigned int min_observation_num_;

  /// Maximum number of observations
  unsigned int max_observation_num_;

  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// Method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////

  /// Load real data settings
  void LoadRealDataSettings(const tinyxml2::XMLElement* settings);

  /// Load the attributes of the cognitive scores
  void LoadRealCognitiveScores(const tinyxml2::XMLElement* settings);

  /// Load the attributes of the landmarks
  void LoadRealLandmarks(const tinyxml2::XMLElement* settings);

  /// Load simulated data settings
  void LoadSimulatedDataSettings(const tinyxml2::XMLElement* settings);

  DataSettings(const DataSettings&);
  DataSettings& operator=(const DataSettings&);
};

}
