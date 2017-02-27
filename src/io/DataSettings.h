#ifndef _DataSettings_h
#define _DataSettings_h

typedef double ScalarType;

#include <string>
#include <fstream>
#include <iostream>
#include <src/observations/Observations.h>

#include "tinyxml2.h"

namespace io {

class DataSettings {
public:


    ////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    DataSettings(const char *XMLFile);

    ~DataSettings();

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    std::string GetPathToGroup() const { return m_PathToGroup; }

    std::string GetPathToTimepoints() const { return m_PathToTimepoints; }

    std::string GetPathToCognitiveScores() const { return m_PathToCognitiveScores; }
    
    std::string GetPathToLandmarks() const { return m_PathToLandmarks; }

    unsigned int GetNumberOfSimulatedSubjects() const { return m_NumberOfSubjects; }

    unsigned int
    GetMinimumNumberOfObservations() const { return m_MinimumNumberOfObservations; }

    unsigned int
    GetMaximumNumberOfObservations() const { return m_MaximumNumberOfObservations; }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    bool IsReal() const { return m_RealData; }
    
    bool LandmarkPresence() const { return m_Landmarks; }
    
    bool CognitiveScoresPresence() const { return m_CognitiveScores; }
        

private:

    /// Read real data if true
    bool m_RealData;
    
    /// Presence of landmarks
    bool m_Landmarks;
    
    /// Presence of cognitive scores
    bool m_CognitiveScores;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Attribute(s) for real data :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    /// Path to data
    std::string m_PathToData;

    /// Path to group csv
    std::string m_PathToGroup;

    /// Path to timepoint csv
    std::string m_PathToTimepoints;

    /// Path to cognitive scores (csv file)
    std::string m_PathToCognitiveScores;
    
    /// Path to landmarks (csv file)
    std::string m_PathToLandmarks;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Attribute(s) for simulated data :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    /// NumberOfIndividuals to simulate
    unsigned int m_NumberOfSubjects;

    /// Minimum number of observations
    unsigned int m_MinimumNumberOfObservations;

    /// Maximum number of observations
    unsigned int m_MaximumNumberOfObservations;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    /// Load real data settings
    void LoadRealDataSettings(const tinyxml2::XMLElement* Settings);

    /// Load simulated data settings
    void LoadSimulatedDataSettings(const tinyxml2::XMLElement* Settings);
    

};

}

#endif //_DataSettings_h
