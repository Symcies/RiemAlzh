#include "AlgorithmSettings.h"

namespace io {


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


AlgorithmSettings::AlgorithmSettings(std::string xml_file) {
  InputsAssert::IsValidAlgoXML(xml_file);

  tinyxml2::XMLDocument parameters;
  parameters.LoadFile(xml_file.c_str());

  auto settings = parameters.FirstChildElement("algorithm-settings");

  max_num_iter_ = std::stoi(settings->FirstChildElement("max-iterations")->GetText());
  num_burnin_iter_ = std::stoi(settings->FirstChildElement("burn-in")->GetText());
  output_iter_ = std::stoi(settings->FirstChildElement("step-size-to-display")->GetText());
  data_save_iter_ = std::stoi(settings->FirstChildElement("step-size-to-save")->GetText());

}


AlgorithmSettings::~AlgorithmSettings() {

}

////////////////////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////////////////////



}
