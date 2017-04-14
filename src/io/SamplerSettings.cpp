#include "SamplerSettings.h"

namespace io {

SamplerSettings::SamplerSettings(const char *xml_file) {
  InputsAssert::IsValidSamplerXML(xml_file);

  number_of_samplers_ = 0;
  sampler_types_.clear();
  sampler_number_of_iterations_.clear();
  
  tinyxml2::XMLDocument parameters;
  parameters.LoadFile(xml_file);

  auto settings = parameters.FirstChildElement("sampler-settings");
  
  for(auto child = settings->FirstChildElement(); child != NULL; child = child->NextSiblingElement()) {
    AddSampler(child);
  }
  
}

SamplerSettings::~SamplerSettings() {}
  


void SamplerSettings::AddSampler(tinyxml2::XMLElement* sampler) {
  std::string type = sampler->FirstChildElement("type")->GetText();
  unsigned int number_of_iter = std::stoul(sampler->FirstChildElement("number-of-iterations")->GetText());
  
  ++number_of_samplers_;
  sampler_types_.push_back(type);
  sampler_number_of_iterations_.push_back(number_of_iter);
}

}