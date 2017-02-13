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
    m_RealData = false;
    
    std::cout << "TODO ... to complete a data class ..." << std::endl;
}
