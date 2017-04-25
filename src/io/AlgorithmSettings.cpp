#include <src/global.h>
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

  int max_num_iter, num_burnin_iter, output_iter, data_save_iter;
  max_num_iter = std::stoi(settings->FirstChildElement("max-iterations")->GetText());
  num_burnin_iter = std::stoi(settings->FirstChildElement("burn-in")->GetText());
  output_iter = std::stoi(settings->FirstChildElement("step-size-to-display")->GetText());
  data_save_iter = std::stoi(settings->FirstChildElement("step-size-to-save")->GetText());

  if(max_num_iter < 0 || num_burnin_iter < 0 || output_iter < 0 || data_save_iter < 0){
    throw InputException("All iterations values must be superior or equal to 0.");
  }

  max_num_iter_ = max_num_iter;
  num_burnin_iter_ = num_burnin_iter;
  output_iter_ = output_iter;
  data_save_iter_ = data_save_iter;
  GV::MAX_ITER = max_num_iter;
}


AlgorithmSettings::~AlgorithmSettings() {

}

////////////////////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////////////////////



}
