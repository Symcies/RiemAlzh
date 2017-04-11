#include "ModelSettings.h"

namespace io {

////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////

ModelSettings::ModelSettings(std::string xml_file) {
  if (InputsAssert::IsValidModelXML(xml_file)){

    tinyxml2::XMLDocument parameters;
    parameters.LoadFile(xml_file.c_str());

    auto settings = parameters.FirstChildElement("model-settings");

    type_ = settings->FirstChildElement("type")->GetText();
    independent_sources_nb_ = std::stoi(settings->FirstChildElement("number-of-independent-sources")->GetText());

    LoadModel(settings);

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

void ModelSettings::LoadModel(const tinyxml2::XMLElement *settings){
  if (InputsAssert::ToLowerCase(type_) == "fastnetwork") {
    LoadFastNetwork(settings);
  } else if (InputsAssert::ToLowerCase(type_) == "meshwork") {
    LoadMeshworkModel(settings);
  } else if (InputsAssert::ToLowerCase(type_) == "network") {
    LoadNetworkModel(settings);
  } else if (InputsAssert::ToLowerCase(type_) == "multivariate") {
    LoadMultivariate(settings);
  } else if (InputsAssert::ToLowerCase(type_) == "univariate") {
    LoadUnivariate(settings);
  } else {
    std::cerr << "The model type defined in model_settings.xml should be one of the following: "
            "Univariate, Multivariate, FastNetwork, Meshwork. Here, it was " << type_ << std::endl ;
  }
}

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
