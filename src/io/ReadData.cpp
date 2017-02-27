#include "ReadData.h"

namespace io {

ReadData::Data
ReadData
::OpenFilesMultivariate(DataSettings &DS, int NbMaxOfSubjects = -1) {
    Data D;

    std::ifstream IndivID(DS.GetPathToGroup());
    std::ifstream DataX(DS.GetPathToTimepoints());
    std::ifstream DataY(DS.GetPathToLandmarks());

    /// Open the Group file;
    if (IndivID.is_open()) {
        int i = 0;
        int Count = 0;
        std::string line;
        std::vector<std::pair<LinearAlgebra<ScalarType>::VectorType, double> > IndivData;
        std::pair<LinearAlgebra<ScalarType>::VectorType, double> Observations;

        while (getline(IndivID, line)) {
            int j = std::stoi(line);
            if (j == i) {
                IndivData.push_back(Observations);
            } else {
                D.push_back(IndivData);
                IndivData.clear();
                IndivData.push_back(Observations);
                i = j;
            }
            if (Count >= NbMaxOfSubjects && NbMaxOfSubjects != -1) {
                break;
            }

        }
        D.push_back(IndivData);
        ++Count;
    } else { std::cout << "Unable to open indiv id's"; }

    /// Open the DATA : timepoints
    if (DataX.is_open()) {
        std::string line;
        getline(DataX, line);
        for (auto it = D.begin(); it != D.end(); it++) {
            for (auto it2 = it->begin(); it2 != it->end(); ++it2) {
                it2->second = std::stod(line);
                getline(DataX, line);
            }
        }
    } else { std::cout << "Unable to open timepoints"; }

    /// Open the DATA : observations
    if (DataY.is_open()) {
        std::string line;
        getline(DataY, line);
        for (auto it = D.begin(); it != D.end(); ++it) {
            for (auto it2 = it->begin(); it2 != it->end(); ++it2) {
                LinearAlgebra<ScalarType>::VectorType X(1827);
                // TODO : check if everything was parsed
                int i = 0;
                std::stringstream LineStream(line);
                std::string cell;
                while (std::getline(LineStream, cell, ',')) {
                    X(i) = std::stod(cell);
                    i++;
                }
                it2->first = X;
                getline(DataY, line);
            }
        }
    } else { std::cout << "Unable to open the observation values"; }

    if (D[0].size() == 0) {
        D.erase(D.begin());
    }

    return D;
}
    
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
        
        /// Same subject
        if(NewSubjectID == CurrentSubjectID)
        {
            TimePoints.push_back(stod(TimePointsLine));
            if(DS.LandmarkPresence()) 
            {
                VectorType NewObs;
                std::stringstream LineStream(LandmarksLine);
                std::string cell;
                while(std::getline(LineStream, cell, ','))
                {
                    NewObs.push_back(std::stod(cell));
                }
                Landmarks.push_back(NewObs);
            }
            if(DS.CognitiveScoresPresence())
            {
                VectorType NewObs;
                std::stringstream LineStream(CognitiveScoresLine);
                std::string cell;
                while(std::getline(LineStream, cell, ','))
                {
                    NewObs.push_back(std::stod(cell));
                }
                CognitiveScores.push_back(NewObs);
            }
        }
        /// New Subject    
        else
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