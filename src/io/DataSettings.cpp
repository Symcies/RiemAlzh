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
  const tinyxml2::XMLElement* Data = Settings->FirstChildElement("observations");
  const tinyxml2::XMLElement* CognitiveScores = Data->FirstChildElement("cognitive-scores");
  const tinyxml2::XMLElement* Landmarks       = Data->FirstChildElement("landmarks");
  
  LoadRealCognitiveScores(CognitiveScores);
  LoadRealLandmarks(Landmarks);

  
}


void
DataSettings
::LoadRealCognitiveScores(const tinyxml2::XMLElement *Settings) {
  
  std::string Presence = Settings->FirstChildElement("presence")->GetText();
  
  if(Presence == "yes") 
  {
    m_CognitiveScoresPresence = true;
    m_PathToCognitiveScores = m_PathToData + Settings->FirstChildElement("path-to-data")->GetText();
    
    if (!std::ifstream(m_PathToCognitiveScores))
    {
      std::cerr << "The file located at " << m_PathToCognitiveScores << " does not exist" << std::endl;
    } 

    m_CognitiveScoresDimension = atoi(Settings->FirstChildElement("dimension")->GetText());
    

    std::cout << "The model is reading cognitive scores located at " << m_PathToCognitiveScores << std::endl;
  }
  else
  {
    m_CognitiveScoresPresence = false;
  }
}

void 
DataSettings
::LoadRealLandmarks(const tinyxml2::XMLElement *Settings) 
{
  std::string Presence = Settings->FirstChildElement("presence")->GetText();
  
  if(Presence == "yes")
  {
    m_LandmarksPresence = true;
    m_PathToLandmarks = m_PathToData + Settings->FirstChildElement("path-to-data")->GetText();
    
    if(!std::ifstream(m_PathToLandmarks))
    {
      std::cerr << "The file located at " << m_PathToLandmarks << " does not exist" << std::endl;
    }
    
    m_LandmarksDimension = atoi(Settings->FirstChildElement("dimension")->GetText());
    
    std::cout << "The model is reading landmarks located at " << m_PathToLandmarks << std::endl;
  }
  else
  {
    m_LandmarksPresence = false;
  }
}


void
DataSettings
::LoadSimulatedDataSettings(const tinyxml2::XMLElement* Settings) {

  m_RealData = false;

  m_NumberOfSubjects = atoi(Settings->FirstChildElement("number-of-individuals")->GetText());
  m_MinimumNumberOfObservations = atoi(Settings->FirstChildElement("min-number-of-observations")->GetText());
  m_MaximumNumberOfObservations = atoi(Settings->FirstChildElement("max-number-of-observations")->GetText());

  std::cout << "The model is simulating between ";
  std::cout << m_MinimumNumberOfObservations << " and " << m_MaximumNumberOfObservations;
  std::cout << " observations for " << m_NumberOfSubjects << " subjects" << std::endl;
}

  
  
}