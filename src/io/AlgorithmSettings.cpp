#include "AlgorithmSettings.h"

namespace io {


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


AlgorithmSettings::AlgorithmSettings(char *xml_file) {
  if (InputsAssert::IsFileCorrect(xml_file, true)){

    tinyxml2::XMLDocument parameters;
    parameters.LoadFile(xml_file);

    auto settings = parameters.FirstChildElement("algorithm-settings");

    max_num_iter_ = std::stoi(settings->FirstChildElement("max-iterations")->GetText());
    num_burnin_iter_ = std::stoi(settings->FirstChildElement("burn-in")->GetText());
    output_iter_ = std::stoi(settings->FirstChildElement("step-size-to-display")->GetText());
    data_save_iter_ = std::stoi(settings->FirstChildElement("step-size-to-save")->GetText());
  }
}


AlgorithmSettings::~AlgorithmSettings() {

}

////////////////////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////////////////////



}
