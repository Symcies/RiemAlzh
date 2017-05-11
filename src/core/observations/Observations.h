#pragma once

typedef double ScalarType;

#include "IndividualObservations.h"
#include "LinearAlgebra.h"

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
  Observations(const Observations &);

  Observations& operator=(const Observations &);

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  unsigned int GetNumberOfSubjects() const { return tot_indiv_num_; }

  ScalarType GetTotalNumberOfObservations() const { return tot_obs_num_; }

  ScalarType GetTotalSumOfCognitiveScores() const { return cog_scores_tot_sum_; }

  ScalarType GetTotalSumOfLandmarks() const { return landmarks_tot_sum_; }

  std::vector<VectorType> GetObservations() const { return indiv_obs_; }
  
  VectorType GetSubjectTimePoints(unsigned int indiv_num) const { return indiv_obs_.at(indiv_num); }

  unsigned int GetNumberOfTimePoints(unsigned int indiv_num) const { return data_.at(indiv_num).GetNumberOfTimePoints(); }

  ScalarType GetSubjectTimePoint(unsigned int indiv_num, unsigned int time_point_num) const {
      return data_.at(indiv_num).GetTimePoint(time_point_num);
  }

  const VectorType& GetSubjectLandmark(unsigned int indiv_num, unsigned int time_point_num) const {
      return data_.at(indiv_num).GetLandmark(time_point_num);
  }

  const VectorType& GetSubjectCognitiveScore(unsigned int indiv_num, unsigned int time_point_num) const {
      return data_.at(indiv_num).GetCognitiveScore(time_point_num);
  }

  const IndividualObservations& GetSubjectObservations(unsigned int indiv_num) const { return data_.at(indiv_num); }

  const int GetId(unsigned int indiv_num) const { return data_.at(indiv_num).GetId(); }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Add the observations of a new patient
  void AddIndividualData(IndividualObservations& indiv_obs);

  /// Add the observations of several new patients
  void AddObservations(Observations& obs);

  /// Initialize the cross-subjects attributes
  void InitializeGlobalAttributes();

private:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Subject attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Individuals
  std::vector<IndividualObservations> data_;

  /// Individual time points
  std::vector<VectorType> indiv_obs_;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Cross-subject attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Total number of subjects
  unsigned int tot_indiv_num_;

  /// Total number of observations
  ScalarType tot_obs_num_;

  /// Sum of the cognitive scores
  ScalarType cog_scores_tot_sum_;

  /// Sum of the landmarks
  ScalarType landmarks_tot_sum_;




};
