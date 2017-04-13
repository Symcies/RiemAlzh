#include "ReadData.h"

namespace io {

Observations ReadData::ReadObservations(const RealDataSettings &ds)
{
  Observations obs;

  InputsAssert::IsFileCorrect(ds.GetPathToGroup(), false);
  InputsAssert::IsFileCorrect(ds.GetPathToTimepoints(), false);
  if (ds.LandmarkPresence()){
    InputsAssert::IsFileCorrect(ds.GetPathToLandmarks(), false);
  }
  if (ds.CognitiveScoresPresence()){
    InputsAssert::IsFileCorrect(ds.GetPathToCognitiveScores(), false);
  }

  std::ifstream indiv_id_file(ds.GetPathToGroup());
  std::ifstream time_points_file(ds.GetPathToTimepoints());
  std::ifstream landmarks_file(ds.GetPathToLandmarks());
  std::ifstream cognitive_scores_file(ds.GetPathToCognitiveScores());

  std::string group_line, time_points_line;
  unsigned int current_subject_id = -1;

  VectorType time_points;
  std::vector<VectorType> landmarks, cognitive_scores;

  while(std::getline(indiv_id_file, group_line))
  {
    if(current_subject_id == -1) {
      current_subject_id = std::stoi(group_line);
    }

    unsigned int new_subject_id = std::stoi(group_line);

    /// If we changed subject
    if(new_subject_id != current_subject_id) {
      IndividualObservations individual = CreateIndividualObs(ds, time_points, landmarks, cognitive_scores);
      obs.AddIndividualData(individual);

      current_subject_id = new_subject_id;
      time_points.clear();
      landmarks.clear();
      cognitive_scores.clear();
    }

    std::getline(time_points_file, time_points_line);
    time_points.push_back(stod(time_points_line));

    if(ds.LandmarkPresence()) {
      VectorType new_obs = ExtractObservation(landmarks_file, ds.GetLandmarksDimension());
      landmarks.push_back(new_obs);
    }

    if(ds.CognitiveScoresPresence()) {
      VectorType new_obs = ExtractObservation(cognitive_scores_file, ds.GetCognitiveScoresDimension());
      cognitive_scores.push_back(new_obs);
    }
  } //end while
  IndividualObservations individual = CreateIndividualObs(ds, time_points, landmarks, cognitive_scores);
  obs.AddIndividualData(individual);

  return obs;
}

IndividualObservations ReadData::CreateIndividualObs(const RealDataSettings& ds, ReadData::VectorType& time_points,
std::vector<ReadData::VectorType>& landmarks, std::vector<ReadData::VectorType>& cognitive_scores){
    IndividualObservations individual(time_points);
    if(ds.LandmarkPresence()) {
      individual.AddLandmarks(landmarks);
    }
    if(ds.CognitiveScoresPresence()) {
      individual.AddCognitiveScores(cognitive_scores);
    }
    return individual;
}

ReadData::VectorType ReadData::ExtractObservation(std::ifstream& file_stream, unsigned int dimension){
  std::string line;
  std::getline(file_stream, line);
  VectorType new_obs(dimension);
  int i = 0;
  std::stringstream line_stream(line);
  std::string cell;
  while(std::getline(line_stream, cell, ',')){
    new_obs(i) = std::stod(cell);
    ++i;
  }
  return new_obs;
}

ReadData::MatrixType ReadData::OpenKernel(std::string file_path) {
  InputsAssert::IsFileCorrect(file_path, false);
  std::ifstream file(file_path);

  //Filling the kernel
  MatrixType kernel(CountLineNumber(file), 1827, 0);

  if (file.is_open()) {
    std::string line;
    int i = 0;
    while (std::getline(file, line)) {
      int j = 0;
      std::stringstream line_stream(line);
      std::string cell;
      while (std::getline(line_stream, cell, ',')) {
        auto a = std::stold(cell);
        if (std::fabs(a) < 10e-13) {
          kernel(i, j) = 0;
        } else {
          kernel(i, j) = std::stold(cell);
        }
        ++j;
      }
      ++i;
    }
  } else {
    std::cout << "Unable to open the kernel located in " + file_path;
  }

  return kernel;
}

int ReadData::CountLineNumber(std::ifstream& file){
  int lines_nb = 0;
  std::string line;
  while (std::getline(file, line)) {
    ++lines_nb;
  }
  file.clear();
  file.seekg(00, std::ios::beg);
  return lines_nb;
}

} //end namespace
