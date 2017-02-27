#include "ModelSettings.h"

namespace io {

////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////

ModelSettings
::ModelSettings(const char *XMLFile) {
    tinyxml2::XMLDocument Parameters;
    Parameters.LoadFile(XMLFile);

    auto Settings = Parameters.FirstChildElement("model-settings");

    std::string Type = Settings->FirstChildElement("type")->GetText();
    m_Type = Type;

    if (Type == "FastNetwork")
        LoadFastNetwork(XMLFile);
    else if (Type == "Meshwork")
        LoadMeshworkModel(XMLFile);
    else if (Type == "Network")
        LoadNetworkModel(XMLFile);
    else
        std::cerr
                << "The model type defined in model_settings.xml should be in {FastNetwork, Meshwork}";
}

ModelSettings
::~ModelSettings() {

}


////////////////////////////////////////////////////////////////////////////////////////////////
/// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////

void
ModelSettings
::LoadFastNetwork(const char *XMLFile) {
    tinyxml2::XMLDocument Parameters;
    Parameters.LoadFile(XMLFile);

    auto Settings = Parameters.FirstChildElement("model-settings");

    m_ManifoldDimension = atoi(Settings->FirstChildElement("number-of-dimensions")->GetText());
    m_NbIndependentSources = atoi(
            Settings->FirstChildElement("number-of-independent-sources")->GetText());

    m_InvertKernelMatrixPath = Settings->FirstChildElement("path-to-kernel-invKd")->GetText();
    m_InterpolationMatrixPath = Settings->FirstChildElement("path-to-kernel-Kxd")->GetText();

    std::cout << "The model used is the " << m_Type << " model." << std::endl;
}


void
ModelSettings
::LoadMeshworkModel(const char *XMLFile) {
    tinyxml2::XMLDocument Parameters;
    Parameters.LoadFile(XMLFile);

    auto Settings = Parameters.FirstChildElement("model-settings");

    m_ManifoldDimension = atoi(Settings->FirstChildElement("number-of-dimensions")->GetText());
    m_NbIndependentSources = atoi(
            Settings->FirstChildElement("number-of-independent-sources")->GetText());

    m_InvertKernelMatrixPath = Settings->FirstChildElement("path-to-kernel-invKd")->GetText();
    m_InterpolationMatrixPath = Settings->FirstChildElement("path-to-kernel-Kxd")->GetText();

    std::cout << "The model used is the " << m_Type << " model." << std::endl;
}

void
ModelSettings
::LoadNetworkModel(const char *XMLFile) {
    tinyxml2::XMLDocument Parameters;
    Parameters.LoadFile(XMLFile);

    auto Settings = Parameters.FirstChildElement("model-settings");

    m_ManifoldDimension = atoi(Settings->FirstChildElement("number-of-dimensions")->GetText());
    m_NbIndependentSources = atoi(
            Settings->FirstChildElement("number-of-independent-sources")->GetText());

    m_InvertKernelMatrixPath = Settings->FirstChildElement("path-to-kernel-invKd")->GetText();
    m_InterpolationMatrixPath = Settings->FirstChildElement("path-to-kernel-Kxd")->GetText();

    std::cout << "The model used is the " << m_Type << " model." << std::endl;
}

}