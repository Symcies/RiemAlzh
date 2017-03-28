#include "ReadData.h"

namespace io {

//TODO(clem): divide and simplify
Observations ReadData::ReadObservations(DataSettings &ds)
{
  Observations obs;

  std::ifstream indiv_id(ds.GetPathToGroup());
  std::ifstream time_points_file(ds.GetPathToTimepoints());
  std::ifstream landmarks_file(ds.GetPathToLandmarks());
  std::ifstream cognitive_scores_file(ds.GetPathToCognitiveScores());

  std::string group_line, time_points_line,
              landmarks_line, cognitive_scores_line;
  unsigned int current_suject_id = -1;

  VectorType time_points;
  std::vector<VectorType> landmarks;
  std::vector<VectorType> cognitive_scores;

  while(std::getline(indiv_id, group_line))
  {
    //current_suject_id initialization
    if(current_suject_id == -1) {
      current_suject_id = std::stoi(group_line);
    }


    unsigned int new_subject_id = std::stoi(group_line);

    /// If we changed subject
    if(new_subject_id != current_suject_id) {
      IndividualObservations individual(time_points);
      if(ds.LandmarkPresence())        individual.AddLandmarks(landmarks);
      if(ds.CognitiveScoresPresence()) individual.AddCognitiveScores(cognitive_scores);
      obs.AddIndividualData(individual);

      current_suject_id = new_subject_id;
      time_points.clear();
      landmarks.clear();
      cognitive_scores.clear();
    }

    std::getline(time_points_file, time_points_line);
    time_points.push_back(stod(time_points_line));

    if(ds.LandmarkPresence()) {
      std::getline(landmarks_file, landmarks_line);
      VectorType new_obs(ds.GetLandmarksDimension());
      int i = 0;
      std::stringstream line_stream(landmarks_line);
      std::string cell;
      while(std::getline(line_stream, cell, ',')) {
        new_obs(i) = std::stod(cell);
        ++i;
      }
      landmarks.push_back(new_obs);
    }

    if(ds.CognitiveScoresPresence()) {
      std::getline(cognitive_scores_file, cognitive_scores_line);
      VectorType new_obs(ds.GetCognitiveScoresDimension());
      int i = 0;
      std::stringstream line_stream(cognitive_scores_line);
      std::string cell;
      while(std::getline(line_stream, cell, ',')){
        new_obs(i) = std::stod(cell);
        ++i;
      }
      cognitive_scores.push_back(new_obs);
    }
  } //end while

  IndividualObservations individual(time_points);
  if(ds.LandmarkPresence())        individual.AddLandmarks(landmarks);
  if(ds.CognitiveScoresPresence()) individual.AddCognitiveScores(cognitive_scores);
  obs.AddIndividualData(individual);

  return obs;
}


ReadData::MatrixType ReadData::OpenKernel(std::string file_path) {
  if(file_path == ""){
    std::cerr << "Error in OpenKernel. Empty file path" << std::endl;
    //TODO: change error return
    return MatrixType(0,0);
  }
  std::ifstream file(file_path);

  //Counting the total number of lines in the file
  int lines_nb = 0;
  std::string line;
  while (std::getline(file, line)) {
    ++lines_nb;
  }
  file.clear();
  file.seekg(00, std::ios::beg);

  //Filling the kernel
  MatrixType kernel(lines_nb, 1827, 0);

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

} //end namespace
