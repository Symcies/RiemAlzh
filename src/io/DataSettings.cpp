#include "DataSettings.h"

namespace io {


////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////

DataSettings
::DataSettings(const char *XMLFile) {
    tinyxml2::XMLDocument File;
    File.LoadFile(XMLFile);

    auto Settings = File.FirstChildElement("data-settings");

    std::string RealData = Settings->FirstChildElement("data")->GetText();

    if (RealData == "true")  LoadRealDataSettings(Settings->FirstChildElement("real-data"));
    if (RealData == "false") LoadSimulatedDataSettings(Settings->FirstChildElement("simulated-data"));
    
}

DataSettings
::~DataSettings() {

}


////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////

void
DataSettings
::LoadRealDataSettings(const tinyxml2::XMLElement* Settings) {
    
    m_RealData = true;
    m_PathToData = Settings->FirstChildElement("folder-path")->GetText();

    m_PathToGroup = m_PathToData + Settings->FirstChildElement("group-file")->GetText();
    m_PathToTimepoints = m_PathToData + Settings->FirstChildElement("timepoints-file")->GetText();

    // TODO : For the path to data, path to timepoints and path to group, need to check if it is ok!!! 
    /// Read the Cognitive scores 
    std::string CognitiveScores = Settings->FirstChildElement("observations")->FirstChildElement("cognitive-scores")->GetText();
    if(CognitiveScores == "none")
    {
        m_CognitiveScores = false;
    }
    else if(!std::ifstream(m_PathToData + CognitiveScores.c_str())) 
    { 
        std::cerr << " The file located at " << CognitiveScores << " does not exist" ; 
    }
    else
    {
        m_CognitiveScores = true;
        m_PathToCognitiveScores = m_PathToData + CognitiveScores;
    }
    
    /// Read the landmarks
    std::string Landmarks = Settings->FirstChildElement("observations")->FirstChildElement("landmarks")->GetText();
    if(Landmarks == "none")
    {
        m_Landmarks = false;
    }
    else if(!std::ifstream(m_PathToData + Landmarks.c_str())) 
    { 
        std::cout << " The file located at " << Landmarks << " does not exist" ; 
    }
    else
    {
        m_Landmarks = true;
        m_PathToLandmarks = m_PathToData + Landmarks;
    }
    

    std::cout << "The model is reading real data located at " << m_PathToData << std::endl;
}


void
DataSettings
::LoadSimulatedDataSettings(const tinyxml2::XMLElement* Settings) {

    m_RealData = false;

    m_NumberOfSubjects = atoi(Settings->FirstChildElement("number-of-individuals")->GetText());
    m_MinimumNumberOfObservations = atoi(
            Settings->FirstChildElement("min-number-of-observations")->GetText());
    m_MaximumNumberOfObservations = atoi(
            Settings->FirstChildElement("max-number-of-observations")->GetText());

    std::cout << "The model is simulating between ";
    std::cout << m_MinimumNumberOfObservations << " and " << m_MaximumNumberOfObservations;
    std::cout << " observations for " << m_NumberOfSubjects << " subjects" << std::endl;
}

    
    
}