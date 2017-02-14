#include "DataSettings.h"


////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////

DataSettings
::DataSettings(const char *XMLFile) 
{
    tinyxml2::XMLDocument File;
    File.LoadFile(XMLFile);
    
    auto Settings = File.FirstChildElement("data-settings");
    
    std::string RealData = Settings->FirstChildElement("real-data")->GetText(); 
    
    if(RealData == "true")  LoadRealDataSettings(XMLFile);
    if(RealData == "false") LoadSimulatedDataSettings(XMLFile);
}

DataSettings
::~DataSettings() 
{
    
}


////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////

void
DataSettings
::LoadRealDataSettings(const char *XMLFile) 
{
    tinyxml2::XMLDocument File;
    File.LoadFile(XMLFile);
    
    auto Settings = File.FirstChildElement("data-settings");
    
    m_RealData = true;
    m_PathToData = Settings->FirstChildElement("folder-path")->GetText(); 
    
    m_PathToGroup        = m_PathToData + Settings->FirstChildElement("group-file")->GetText();
    m_PathToTimepoints   = m_PathToData + Settings->FirstChildElement("timepoints-file")->GetText();
    m_PathToObservations = m_PathToData + Settings->FirstChildElement("observations-file")->GetText();
    
    std::cout << "The model is reading real data located at " << m_PathToData << std::endl;
}


void
DataSettings
::LoadSimulatedDataSettings(const char *XMLFile) 
{
    tinyxml2::XMLDocument File;
    File.LoadFile(XMLFile);
    
    auto Settings = File.FirstChildElement("data-settings");
    
    m_RealData = false;
    
    m_NumberOfSubjects = atoi(Settings->FirstChildElement("number-of-individuals")->GetText());
    m_MinimumNumberOfObservations = atoi(Settings->FirstChildElement("min-number-of-observations")->GetText());
    m_MaximumNumberOfObservations = atoi(Settings->FirstChildElement("max-number-of-observations")->GetText());
    
    std::cout << "The model is simulating the between ";
    std::cout << m_MinimumNumberOfObservations << " and " << m_MaximumNumberOfObservations;
    std::cout << " for " << m_NumberOfSubjects << " subjects" << std::endl;
}
