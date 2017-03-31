#include "ModelSettings.h"

namespace io {

////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////

ModelSettings::ModelSettings(const char *xml_file) {
  tinyxml2::XMLDocument parameters;
  parameters.LoadFile(xml_file);

  auto settings = parameters.FirstChildElement("model-settings");

  type_ = settings->FirstChildElement("type")->GetText();
  independent_sources_nb_ = atoi(settings->FirstChildElement("number-of-independent-sources")->GetText());

  if (type_ == "FastNetwork") {
    LoadFastNetwork(settings);
  } else if (type_ == "Meshwork") {
    LoadMeshworkModel(settings);
  } else if (type_ == "Network") {
    LoadNetworkModel(settings);
  } else if (type_ == "Multivariate") {
    LoadMultivariate(settings);
  }else if (type_ == "Univariate") {
    LoadUnivariate(settings);
  } else {
    std::cerr << "The model type defined in model_settings.xml should be in {FastNetwork, Meshwork}";
  }
}

ModelSettings::~ModelSettings() {
}


////////////////////////////////////////////////////////////////////////////////////////////////
/// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////

void ModelSettings::LoadFastNetwork(const tinyxml2::XMLElement *settings) {
  invert_kernel_matrix_path_ = settings->FirstChildElement("path-to-kernel-invKd")->GetText();
  interpolation_matrix_path_ = settings->FirstChildElement("path-to-kernel-Kxd")->GetText();

  PrintModelInfo();
}


void ModelSettings::LoadMeshworkModel(const tinyxml2::XMLElement *settings) {
  invert_kernel_matrix_path_ = settings->FirstChildElement("path-to-kernel-invKd")->GetText();
  interpolation_matrix_path_ = settings->FirstChildElement("path-to-kernel-Kxd")->GetText();

  PrintModelInfo();
}

void ModelSettings::LoadNetworkModel(const tinyxml2::XMLElement *settings) {
  invert_kernel_matrix_path_ = settings->FirstChildElement("path-to-kernel-invKd")->GetText();
  interpolation_matrix_path_ = settings->FirstChildElement("path-to-kernel-Kxd")->GetText();

  PrintModelInfo();
}

void ModelSettings::LoadMultivariate(const tinyxml2::XMLElement *settings) {
  PrintModelInfo();
}

void ModelSettings::LoadUnivariate(const tinyxml2::XMLElement *settings) {
  PrintModelInfo();
}

void ModelSettings::PrintModelInfo(){
  std::cout << "The model used is the " << type_ << " model." << std::endl;
}

}
