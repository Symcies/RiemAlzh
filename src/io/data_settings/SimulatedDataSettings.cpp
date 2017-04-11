#include "SimulatedDataSettings.h"


namespace  io {


////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////

SimulatedDataSettings::SimulatedDataSettings(std::string xml_file) : DataSettings(xml_file) {

  InputsAssert::IsValidDataXML(xml_file);

  tinyxml2::XMLDocument file;
  file.LoadFile(xml_file.c_str());
  auto settings = file.FirstChildElement("data-settings")->FirstChildElement("simulated-data");

  subjects_total_num_ = std::stoi(settings->FirstChildElement("number-of-individuals")->GetText());
  min_observation_num_ = std::stoi(settings->FirstChildElement("min-number-of-observations")->GetText());
  max_observation_num_ = std::stoi(settings->FirstChildElement("max-number-of-observations")->GetText());
  dimension_ = std::stoi(settings->FirstChildElement("dimension")->GetText());


  std::cout << "The model is simulating between ";
  std::cout << min_observation_num_ << " and " << max_observation_num_;
  std::cout << " observations for " << subjects_total_num_ << " subjects" << std::endl;


}

  

}