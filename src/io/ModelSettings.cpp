#include "ModelSettings.h"

namespace io {

////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////

ModelSettings::ModelSettings(std::string xml_file) {
  InputsAssert::IsValidModelXML(xml_file);

  tinyxml2::XMLDocument parameters;
  parameters.LoadFile(xml_file.c_str());

  auto settings = parameters.FirstChildElement("model-settings");

  type_ = settings->FirstChildElement("type")->GetText();
  independent_sources_nb_ = std::stoi(settings->FirstChildElement("number-of-independent-sources")->GetText());
  LoadInitialRandomVariables(settings->FirstChildElement("variables"));

  LoadModel(settings);

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

void ModelSettings::LoadInitialRandomVariables(const tinyxml2::XMLElement *settings) {
  for(auto child = settings->FirstChildElement(); child != NULL; child = child->NextSiblingElement()) {
    std::string name = child->FirstChildElement("name")->GetText();
    std::vector<double> initial_params = LoadRVParameters(child->FirstChildElement("initial-parameters"));
    std::vector<double> second_params  = LoadRVParameters(child->FirstChildElement("second-parameters"));
    ScalarType proposition_variance = std::stod(child->FirstChildElement("proposition-variance")->GetText());
    
    auto rv  = std::make_pair(initial_params, proposition_variance);
    auto rv2 = std::make_pair(second_params, proposition_variance);
    
    init_random_variables_[name] = rv;
    second_random_variables_[name] = rv2;
  }
}

std::vector<double> ModelSettings::LoadRVParameters(const tinyxml2::XMLElement* parameters) {
   
  double mean = std::stod(parameters->FirstChildElement("mean")->GetText());
  double variance = std::stod(parameters->FirstChildElement("variance")->GetText());

  return {mean, variance};
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
