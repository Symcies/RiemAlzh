#include "InputsAssert.h"


void InputsAssert::IsFilePathCorrect(std::string path){
  std::ifstream test(path);
  if (test) {
    return;
  }
  throw InputException("The file path '" + path + "' is incorrect.");
}

void InputsAssert::IsFileEmpty(std::string path){
  std::ifstream test(path);
  test.seekg(0, std::ios::end);
  if (test.tellg() == 0){
    throw InputException("The file at '" + path + "' is empty.");
  }
  return;
}

void InputsAssert::IsXMLValid(std::string path){
  tinyxml2::XMLDocument file;
  file.LoadFile(path.c_str());

  if (file.Error()){
    throw InputException("The XML of file " + path + " is not correct.");
  }
}


void InputsAssert::IsValidModelXML(std::string path){
  IsFileCorrect(path, true);

  tinyxml2::XMLDocument file;
  file.LoadFile(path.c_str());

  auto settings = file.FirstChildElement("model-settings");
  if (settings == NULL) {
    throw InputException("The model xml misses the parameter model-settings.");
  }

  auto type = settings->FirstChildElement("type");
  if (type == NULL) {
    throw InputException("The model xml misses the parameter type, "
                           "child of the parameter model-settings.");
  }

  auto indep_source_num = settings->FirstChildElement("number-of-independent-sources");
  if (indep_source_num == NULL) {
    throw InputException("The model xml misses the parameter number-of-independent-sources, "
                                   "child of the parameter model-settings.");
  }

  return;
}

void InputsAssert::IsValidAlgoXML(std::string path){
  IsFileCorrect(path, true);

  tinyxml2::XMLDocument file;
  file.LoadFile(path.c_str());

  auto settings = file.FirstChildElement("algorithm-settings");
  if (settings == NULL) {
    throw InputException("The algorithm xml misses the parameter algorithm-settings.");
  }

  std::string first_order_children[] =
          {"max-iterations", "burn-in", "step-size-to-display", "step-size-to-save"};
  for (int i = 0; i < 4; i++){
    auto child = settings->FirstChildElement(first_order_children[i].c_str());
    if (child == NULL) {
      throw InputException("The algorithm xml misses the parameter " + first_order_children[i] +
                                     ", child of the parameter algorithm-settings.");
    }
  }
};


void InputsAssert::IsValidDataXML(std::string path){
  IsFileCorrect(path, true);

  tinyxml2::XMLDocument file;
  file.LoadFile(path.c_str());

  auto settings = file.FirstChildElement("data-settings");
  if (settings == NULL) {
    throw InputException("The data xml misses the parameter data-settings.");
  }
  auto type = settings->FirstChildElement("data-type");
  if (type == NULL) {
    throw InputException("The data xml misses the parameter data-type, child of data-settings.");
  }

 if (StringToBool(type->GetText())){
    IsValidRealData(settings);
  }
  else {
    IsValidSimulatedData(settings);
  }

};

void InputsAssert::IsValidSimulatedData(const tinyxml2::XMLElement * settings){
  auto simulated = settings->FirstChildElement("simulated-data");
  if (simulated == NULL) {
    throw InputException("The data xml misses the parameter simulated-data, child of data-settings..");
  }

  std::string second_order_children[] =
          {"dimension", "number-of-individuals", "min-number-of-observations", "max-number-of-observations"};
  for (int i = 0; i < 4; i++){
    auto child = simulated->FirstChildElement(second_order_children[i].c_str());
    if (child == NULL) {
      throw InputException("The data xml misses the parameter " + second_order_children[i] +
                             ", child of the parameter simulated-data.");
    }
  }
}

void InputsAssert::IsValidRealData(const tinyxml2::XMLElement * settings){
  auto real = settings->FirstChildElement("real-data");
  if (settings == NULL) {
    throw InputException("The data xml misses the parameter real-data, child of data-settings..");
  }


  std::string second_order_children[] =
          {"folder-path", "group-file", "timepoints-file"};
  for (int i = 0; i < 3; i++){
    auto child = real->FirstChildElement(second_order_children[i].c_str());
    if (child == NULL) {
      throw InputException("The data xml misses the parameter " + second_order_children[i] +
                             ", child of the parameter real-data.");
    }
  }

  auto observations = real->FirstChildElement("observations");
  if (observations == NULL) {
    throw InputException("The data xml misses the parameter observations, child of data-settings..");
  }

  auto cognitive_scores = observations->FirstChildElement("cognitive-scores");
  if (cognitive_scores == NULL) {
    throw InputException("The data xml misses the parameter cognitive-scores, child of observations.");
  }

  std::string fourth_order_children[] =
          {"presence", "path-to-data", "dimension"};
  for (int i = 0; i < 3; i++){
    auto child = cognitive_scores->FirstChildElement(fourth_order_children[i].c_str());
    if (child == NULL) {
      throw InputException("The data xml misses the parameter " + second_order_children[i] +
                             ", child of the parameter cognitive-scores.");
    }
  }

  auto landmarks = observations->FirstChildElement("landmarks");
  if (landmarks == NULL) {
    throw InputException("The data xml misses the parameter landmarks, child of observations.");
  }

  for (int i = 0; i < 3; i++){
    auto child = landmarks->FirstChildElement(fourth_order_children[i].c_str());
    if (child == NULL) {
      throw InputException("The data xml misses the parameter " + second_order_children[i] +
                             ", child of the parameter landmarks.");
    }
  }

}

void InputsAssert::IsFileCorrect(std::string path, bool is_xml){
  IsFilePathCorrect(path);
  IsFileEmpty(path);
  if (is_xml) {
    IsXMLValid(&path[0]);
  }
}

std::string InputsAssert::ToLowerCase(std::string data){
  std::transform(data.begin(), data.end(), data.begin(), ::tolower);
  return data;
}

bool InputsAssert::StringToBool(std::string data){
  std::transform(data.begin(), data.end(), data.begin(), ::tolower);
  if (data == "true") {
    return true;
  } else if (data == "false") {
    return false;
  }
  std::cerr << data << "should be 'true' or 'false'." << std::endl;
}