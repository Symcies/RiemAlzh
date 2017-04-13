#pragma once

#include "DataSettings.h"

namespace io {

class SimulatedDataSettings : public DataSettings {
 
 public:
  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////

  SimulatedDataSettings(std::string xml_file);
  
  
  unsigned int GetDimensionOfSimulatedObservations() const { return dimension_; }
  
  unsigned int GetNumberOfSimulatedSubjects() const { return subjects_total_num_; }

  unsigned int GetMinimumNumberOfObservations() const { return min_observation_num_; }

  unsigned int GetMaximumNumberOfObservations() const { return max_observation_num_; }
  
  
 private:
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// Attribute(s)  :
  ////////////////////////////////////////////////////////////////////////////////////////////////

  /// NumberOfIndividuals to simulate
  unsigned int subjects_total_num_;

  /// Minimum number of observations
  unsigned int min_observation_num_;

  /// Maximum number of observations
  unsigned int max_observation_num_;
  
  /// Dimension of the generated data
  unsigned int dimension_;
  
  
};

}
