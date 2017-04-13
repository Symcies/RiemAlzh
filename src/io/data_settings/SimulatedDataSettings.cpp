#include "SimulatedDataSettings.h"


namespace  io {


////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////

SimulatedDataSettings::SimulatedDataSettings(const char *xml_file) : DataSettings(xml_file) {
  
  tinyxml2::XMLDocument file;
  file.LoadFile(xml_file);
  auto settings = file.FirstChildElement("data-settings")->FirstChildElement("simulated-data");
  

  subjects_total_num_  = atoi(settings->FirstChildElement("number-of-individuals")->GetText());
  min_observation_num_ = atoi(settings->FirstChildElement("min-number-of-observations")->GetText());
  max_observation_num_ = atoi(settings->FirstChildElement("max-number-of-observations")->GetText());
  dimension_ = atoi(settings->FirstChildElement("dimension")->GetText());


  std::cout << "The model is simulating between ";
  std::cout << min_observation_num_ << " and " << max_observation_num_;
  std::cout << " observations for " << subjects_total_num_ << " subjects" << std::endl;


}

  

}