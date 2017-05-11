#include "InputsAssert.h"

/**
 *  @file    InputsAssert.cpp
 *  @date    2017/05/03
 *
 *  @brief Utils file to assert input validity, in terms of path, arguments type, ...
 */

/**
 * @brief Checks is the file path leads to an openable file
 *
 * @param path: the path to check
 * @throw InputException if the file is not openable
 *
 * @TODO check for permissions
 */
void InputsAssert::IsFilePathCorrect(std::string path){
  std::ifstream test(path);
  if (test) {
    return;
  }
  throw InputException("The file path '" + path + "' is incorrect.");
}

/**
 * @brief Checks is the file path leads to an empty file
 *
 * @param path: the path to check
 * @throw InputException if the file is empty
 */
void InputsAssert::IsFileEmpty(std::string path){
  std::ifstream test(path);
  test.seekg(0, std::ios::end);
  if (test.tellg() == 0){
    throw InputException("The file at '" + path + "' is empty.");
  }
  return;
}

/**
 * @brief Checks is the file has a valid XML syntax
 *
 * @param path: the path of the file to check
 * @throw InputException if the XML syntax is wrong
 */
void InputsAssert::IsXMLValid(std::string path){
  IsXMLFilePath(path);

  tinyxml2::XMLDocument file;
  file.LoadFile(path.c_str());

  if (file.Error()){
    throw InputException("The XML of file " + path + " is not correct.");
  }
}

/**
 * @brief Checks is the file has a correct XML for model settings
 *
 * @param path: the path to check
 * @throw InputException if parameters of the XML are missing
 */
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
  CheckElemContentIsString(type);

  auto indep_source_num = settings->FirstChildElement("number-of-independent-sources");
  if (indep_source_num == NULL) {
    throw InputException("The model xml misses the parameter number-of-independent-sources, "
                                 "child of the parameter model-settings.");
  }
  CheckElemContentIsInt(indep_source_num);

  auto variables = settings->FirstChildElement("variables");
  if (variables == NULL) {
    throw InputException("The model xml misses the parameter variables, "
                                 "child of the parameter model-settings.");
  }

  //TODO: to complete
  std::unordered_set<std::string> names;
  for(auto second_order_child = variables->FirstChildElement(); second_order_child != NULL;
      second_order_child = second_order_child->NextSiblingElement()) {
    /// Name
    auto name = second_order_child->FirstChildElement("name");
      if (name == NULL) {
        throw InputException("The model xml misses the parameter name");
      }

    /// Double children
    std::string third_order_children_double[] = {"proposition-variance"};
    names.insert(second_order_child->Value());
    for (int i = 0; i < 1; i++) {
      auto child = second_order_child->FirstChildElement(third_order_children_double[i].c_str());
      if (child == NULL) {
        throw InputException("The model xml misses the parameter " + third_order_children_double[i]);
      }
      CheckElemContentIsDouble(child);
    }

    std::string third_order_children_param[] = {"initial-parameters", "second-parameters"};
    names.insert(second_order_child->Value());
    for (int i = 0; i < 2; i++) {
      auto child_order3 = second_order_child->FirstChildElement(third_order_children_param[i].c_str());
      if (child_order3 == NULL) {
        throw InputException("The model xml misses the parameter " + third_order_children_param[i]);
      }
      std::string fourth_order_children_param[] = {"mean", "variance"};
      for (int i = 0; i < 2; i++) {
        auto child_order4 = child_order3->FirstChildElement(fourth_order_children_param[i].c_str());
        if (child_order4 == NULL) {
          throw InputException("The model xml misses the parameter " + fourth_order_children_param[i]);
        }
        CheckElemContentIsDouble(child_order4);
      }
    }
  }

  AllVariablesPresentInModelXML(type->GetText(), names);

  return;
}

/**
 *
 */
void InputsAssert::AllVariablesPresentInModelXML(std::string type, std::unordered_set<std::string> names){
  if(ToLowerCase(type) == "univariate"){
    std::string univ_names[] = {"noise", "P", "Ksi", "Tau"};
    for(int i = 0; i< 4; i++){
      if (names.count(univ_names[i].c_str())!=1){
        throw InputException("The model xml misses the parameter " + univ_names[i]);
      }
    }

  }
  else if(ToLowerCase(type) == "multivariate"){
    std::string multiv_names[] = {"noise", "G", "Ksi", "Tau", "Delta", "Beta", "S"};
    for(int i = 0; i< 7; i++){
      if (names.count(multiv_names[i].c_str())!=1){
        throw InputException("The model xml misses the parameter " + multiv_names[i]);
      }
    }
  }
}

/**
 * @brief Checks is the file has a correct XML for algorithm settings
 *
 * @param path: the path to check
 * @throw InputException if parameters of the XML are missing
 */
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
    CheckElemContentIsInt(child);
  }
};

/**
 * @brief Checks is the file has a correct XML for sampler settings
 *
 * @param path: the path to check
 * @throw InputException if parameters of the XML are missing
 */
void InputsAssert::IsValidSamplerXML(std::string path){
  IsFileCorrect(path, true);

  tinyxml2::XMLDocument file;
  file.LoadFile(path.c_str());

  auto settings = file.FirstChildElement("sampler-settings");
  if (settings == NULL) {
    throw InputException("The sampler xml misses the parameter sampler-settings.");
  }

  for(auto second_order_child = settings->FirstChildElement(); second_order_child != NULL;
      second_order_child = second_order_child->NextSiblingElement()) {
    /// String children
    std::string third_order_children_str[] = {"type"};
    for (int i = 0; i < 1; i++) {
      auto child = second_order_child->FirstChildElement(third_order_children_str[i].c_str());
      if (child == NULL) {
        throw InputException("The algorithm xml misses the parameter " + third_order_children_str[i] +
                             ", child of the parameter algorithm-settings.");
      }
      CheckElemContentIsString(child);
    }

    /// Int children
    std::string third_order_children_int[] = {"number-of-iterations"};
    for (int i = 0; i < 1; i++) {
      auto child = second_order_child->FirstChildElement(third_order_children_int[i].c_str());
      if (child == NULL) {
        throw InputException("The algorithm xml misses the parameter " + third_order_children_int[i] +
                             ", child of the parameter algorithm-settings.");
      }
      CheckElemContentIsInt(child);
    }
  }

};

/**
 * @brief Checks is the file has a correct XML for data settings
 *
 * @param path: the path to check
 * @throw InputException if parameters of the XML are missing
 */
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
  CheckElemContentIsBool(type);

  if (StringToBool(type->GetText())){
    IsValidRealData(settings);
  }
  else {
    IsValidSimulatedData(settings);
  }

};

/**
 * @brief Checks is the XMLElement has a correct shape for simulated data settings
 *
 * @param settings: the XMLElement to check
 * @throw InputException if parameters of the XML are missing
 */
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

/**
 * @brief Checks is the XMLElement has a correct shape for real data settings
 *
 * @param settings: the XMLElement to check
 * @throw InputException if parameters of the XML are missing
 */
void InputsAssert::IsValidRealData(const tinyxml2::XMLElement * settings){
  auto real = settings->FirstChildElement("real-data");
  if (real == NULL) {
    throw InputException("The data xml misses the parameter real-data, child of data-settings..");
  }

  auto folder = real->FirstChildElement("folder-path");
  if (folder == NULL) {
    throw InputException("The data xml misses the parameter folder-path"
                         ", child of the parameter real-data.");
  }

  std::string second_order_children[] =
          {"group-file", "timepoints-file"};
  for (int i = 0; i < 2; i++){
    auto child = real->FirstChildElement(second_order_children[i].c_str());
    if (child == NULL) {
      throw InputException("The data xml misses the parameter " + second_order_children[i] +
                             ", child of the parameter real-data.");
    }
    IsCsvFilePath(child->GetText());
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

/**
 * @brief Checks is the file is correct (exists, is not empty, and has a correct XML shape)
 *
 * @param path: the path of the file to check
 * @param is_xml: Boolean, which should be true if the file is an xml and false otherwise
 * @throw InputException if the file is incorrect (not existing, empty, or with an incorrect shape)
 */
void InputsAssert::IsFileCorrect(std::string path, bool is_xml){
  IsFilePathCorrect(path);
  IsFileEmpty(path);
  if (is_xml) {
    IsXMLValid(&path[0]);
  }
}

/**
 * @brief Checks is the path has the correct shape and ends with csv (with regex)
 *
 * @param path: the path of the file to check
 * @throw InputException if the path doesn't have the right shape
 */
void InputsAssert::IsCsvFilePath(std::string path){
  if (!std::regex_match(path, std::regex("([\\w\\/\\.]+)\\.csv"))){
    throw InputException("The data xml provides wrong csv file name: " + path);
  }
}

/**
 * @brief Checks is the path has the correct shape and ends with xml (with regex)
 *
 * @param path: the path of the file to check
 * @throw InputException if the path doesn't have the right shape
 */
void InputsAssert::IsXMLFilePath(std::string path){
  if (!std::regex_match(path, std::regex("([\\w\\/\\.]+)\\.xml"))){
    throw InputException("Incorrect path format: " + path);
  }
}

/**
 * @brief Util function to transform a string to lower case
 *
 * @param data: string to convert to lower case
 * @return data string in lower case
 */
std::string InputsAssert::ToLowerCase(std::string data){
  std::transform(data.begin(), data.end(), data.begin(), ::tolower);
  return data;
}

/**
 * @brief Util function to transform a string to a boolean
 *
 * @param data: string to convert to a boolean (can be true, True, t, T or false, False, f, F)
 * @return data in boolean form
 * @throw InputException if the string cannot be converted
 */
bool InputsAssert::StringToBool(std::string data){
  if (ToLowerCase(data) == "true" || ToLowerCase(data) == "t") {
    return true;
  } else if (ToLowerCase(data) == "false" || ToLowerCase(data) == "f") {
    return false;
  }
  throw InputException(data + " should be 'true' or 'false'.");
}

/**
 * @brief Util function which checks if an XMLElement is an integer
 *
 * @param elem: XMLElement to check
 * @throw InputException if the XMLElement is not an integer, or not a leaf
 */
void InputsAssert::CheckElemContentIsInt(const tinyxml2::XMLElement* elem){
  if (elem->FirstChild()->ToText() != NULL) {
    try {
      int res = std::stoi(elem->GetText());
    } catch (std::exception exception) {
      std::string message = std::string("The parameter ") + elem->Value() + std::string(" is not an integer.");
      throw InputException(message);
    }
  }
  else{
    std::string message = std::string("The parameter ") + elem->Value() + std::string(" should be a leaf.");
    throw InputException(message);
  }
}

/**
 * @brief Util function which checks if an XMLElement is a double
 *
 * @param elem: XMLElement to check
 * @throw InputException if the XMLElement is not a double, or not a leaf
 */
void InputsAssert::CheckElemContentIsDouble(const tinyxml2::XMLElement* elem){
  if (elem->FirstChild()->ToText() != NULL) {
    try {
      double res = std::stod(elem->GetText());
    } catch (std::exception exception) {
      std::string message = std::string("The parameter ") + elem->Value() + std::string(" is not a double.");
      throw InputException(message);
    }
  }
  else{
    std::string message = std::string("The parameter ") + elem->Value() + std::string(" should be a leaf.");
    throw InputException(message);
  }
}

/**
 * @brief Util function which checks if an XMLElement is a string
 *
 * @param elem: XMLElement to check
 * @throw InputException if the XMLElement is not a string, or not a leaf
 */
void InputsAssert::CheckElemContentIsString(const tinyxml2::XMLElement* elem){
  if (elem->FirstChild()->ToText() != NULL) {
    try {
      std::string res = elem->GetText();
    } catch (std::exception exception) {
      std::string message = std::string("The parameter ") + elem->Value() + std::string(" is not a string.");
      throw InputException(message);
    }
  }
  else{
    std::string message = std::string("The parameter ") + elem->Value() + std::string(" should be a leaf.");
    throw InputException(message);
  }
}

/**
 * @brief Util function which checks if an XMLElement is a boolean
 *
 * @param elem: XMLElement to check
 * @throw InputException if the XMLElement is not a boolean, or not a leaf
 */
void InputsAssert::CheckElemContentIsBool(const tinyxml2::XMLElement* elem){
  if (elem->FirstChild()->ToText() != NULL) {
    try {
      StringToBool(elem->GetText());
    } catch (std::exception exception) {
      std::string message = std::string("The parameter ") + elem->Value() + std::string(" is not a double.");
      throw InputException(message);
    }
  }
  else{
    std::string message = std::string("The parameter ") + elem->Value() + std::string(" should be a leaf.");
    throw InputException(message);
  }
}