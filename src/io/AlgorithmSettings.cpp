#include "AlgorithmSettings.h"

namespace io {


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


AlgorithmSettings
::AlgorithmSettings(char *XMLFile) {
    tinyxml2::XMLDocument Parameters;
    Parameters.LoadFile(XMLFile);

    auto Settings = Parameters.FirstChildElement("algorithm-settings");

    m_MaximumNumberOfIterations = atoi(
            Settings->FirstChildElement("max-iterations")->GetText());
    m_NumberOfBurnInIterations = atoi(Settings->FirstChildElement("burn-in")->GetText());
    m_CounterToDisplayOutputs = atoi(
            Settings->FirstChildElement("step-size-to-display")->GetText());
    m_CounterToSaveData = atoi(Settings->FirstChildElement("step-size-to-save")->GetText());
}


AlgorithmSettings
::~AlgorithmSettings() {

}

////////////////////////////////////////////////////////////////////////////////////////////////
/// 
////////////////////////////////////////////////////////////////////////////////////////////////



}