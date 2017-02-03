
#include "ModelSettings.h"


    
////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////

ModelSettings
::ModelSettings(const char *XMLFile) 
{
    tinyxml2::XMLDocument Parameters;
    Parameters.LoadFile(XMLFile);
    
    auto Settings = Parameters.FirstChildElement("model-settings");
    
    std::string Type = Settings->FirstChildElement("type")->GetText();
    
    if(Type == "FastNetwork")
        LoadFastNetwork(XMLFile);
    else if(Type == "Meshwork")
        LoadMeshworkModel(XMLFile);
    else 
        std::cerr << "The model type defined in model_settings.xml should be in {FastNetwork, Meshwork}";
}

ModelSettings
::~ModelSettings() 
{
    
}


////////////////////////////////////////////////////////////////////////////////////////////////
/// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////

void
ModelSettings
::LoadFastNetwork(const char *XMLFile) 
{
    
}


void 
ModelSettings
::LoadMeshworkModel(const char *XMLFile) 
{
    tinyxml2::XMLDocument Parameters;
    Parameters.LoadFile(XMLFile);
    
    auto Settings = Parameters.FirstChildElement("model-settings");
    
    m_ManifoldDimension = atoi(Settings->FirstChildElement("manifold-dimension")->GetText());
    m_NbIndependentSources = atoi(Settings->FirstChildElement("number-of-indepenent-sources")->GetText());
    
}