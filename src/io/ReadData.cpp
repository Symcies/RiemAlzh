#include "ReadData.h"

namespace io {

  
Observations
ReadData
::ReadObservations(DataSettings &DS) 
{
  Observations Obs;
  
  std::ifstream IndivID(DS.GetPathToGroup());
  std::ifstream TimePointsFile(DS.GetPathToTimepoints());
  std::ifstream LandmarksFile(DS.GetPathToLandmarks());
  std::ifstream CognitiveScoresFile(DS.GetPathToCognitiveScores());
      
  std::string GroupLine, TimePointsLine, LandmarksLine, CognitiveScoresLine;
  unsigned int CurrentSubjectID = -1;
  
  VectorType TimePoints;
  std::vector<VectorType> Landmarks;
  std::vector<VectorType> CognitiveScores;
  
  while(getline(IndivID, GroupLine))
  {
    if(CurrentSubjectID == -1) CurrentSubjectID = std::stoi(GroupLine);
    
    unsigned int NewSubjectID = std::stoi(GroupLine);
    getline(TimePointsFile, TimePointsLine);
    if(DS.LandmarkPresence())        getline(LandmarksFile, LandmarksLine);
    if(DS.CognitiveScoresPresence()) getline(CognitiveScoresFile, CognitiveScoresLine);
    
    /// New subject
    if(NewSubjectID != CurrentSubjectID)
    {
      IndividualObservations Individual(TimePoints);
      if(DS.LandmarkPresence())        Individual.AddLandmarks(Landmarks);
      if(DS.CognitiveScoresPresence()) Individual.AddCognitiveScores(CognitiveScores);
      Obs.AddIndividualData(Individual);
                  
      CurrentSubjectID = NewSubjectID;
      TimePoints.clear();
      Landmarks.clear();
      CognitiveScores.clear();
    }

    TimePoints.push_back(stod(TimePointsLine));
    if(DS.LandmarkPresence()) 
    {
      VectorType NewObs(DS.GetLandmarksDimension());
      int i = 0;
      std::stringstream LineStream(LandmarksLine);
      std::string cell;
      while(std::getline(LineStream, cell, ','))
      {
        NewObs(i) = std::stod(cell);
        ++i;
      }
      
      Landmarks.push_back(NewObs);
    }
    if(DS.CognitiveScoresPresence())
    {
        
      VectorType NewObs(DS.GetCognitiveScoresDimension());
      int i = 0;
      std::stringstream LineStream(CognitiveScoresLine);
      std::string cell;
      while(std::getline(LineStream, cell, ','))
      {
        NewObs(i) = std::stod(cell);
        ++i;
      }
      CognitiveScores.push_back(NewObs);
    }
  }
  
  IndividualObservations Individual(TimePoints);
  if(DS.LandmarkPresence())        Individual.AddLandmarks(Landmarks);
  if(DS.CognitiveScoresPresence()) Individual.AddCognitiveScores(CognitiveScores);
  Obs.AddIndividualData(Individual);
  
  
  return Obs;
}


ReadData::MatrixType
ReadData
::OpenKernel(std::string FilePath) {
  std::ifstream file(FilePath);
  std::string line;
  int NbLines = 0;
  while (std::getline(file, line)) {
    ++NbLines;
  }
  file.clear();
  file.seekg(00, std::ios::beg);


  MatrixType Kernel(NbLines, 1827, 0);


  if (file.is_open()) {
    std::string line;
    int i = 0;
    while (std::getline(file, line)) {
      int j = 0;
      std::stringstream LineStream(line);
      std::string cell;
      while (std::getline(LineStream, cell, ',')) {
        auto a = std::stold(cell);
        if (fabs(a) < 10e-13) {
          Kernel(i, j) = 0;
        } else {
          Kernel(i, j) = std::stold(cell);
        }

        ++j;
      }
      ++i;
    }
  } else { std::cout << "Unable to open the kernel located in " + FilePath; }


  return Kernel;
}

}