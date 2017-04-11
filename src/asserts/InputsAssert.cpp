#include "InputsAssert.h"


bool InputsAssert::IsFilePathCorrect(std::string path){
  std::ifstream test(path);
  if (test) {
    return true;
  }
  std::cerr << "The file path " << path << " is incorrect." << std::endl;
  throw std::exception();
}

bool InputsAssert::IsFileEmpty(std::string path){
  std::ifstream test(path);
  test.seekg(0, std::ios::end);
  if (test.tellg() == 0){
    std::cerr << "The file " << path << " is empty." << std::endl;
    throw std::exception();
  }
  return false;
}

bool InputsAssert::IsXMLValid(std::string path){
  tinyxml2::XMLDocument file;
  file.LoadFile(path.c_str());

  if (file.Error()){
    throw std::logic_error("The XML of file " + path + " is not correct.");
  }
  return !file.Error();
}


bool InputsAssert::IsValidModelXML(std::string path){
  IsFileCorrect(path, true);

  tinyxml2::XMLDocument file;
  file.LoadFile(path.c_str());

  auto settings = file.FirstChildElement("model-settings");
  if (settings == NULL) {
    std::cerr << "The model xml misses the parameter model-settings." << std::endl;
    throw std::logic_error("Missing parameter in model XML");
  }

  auto type = settings->FirstChildElement("type");
  if (type == NULL) {
    std::cerr << "The model xml misses the parameter type, "
            "child of the parameter model-settings." << std::endl;
    throw std::logic_error("Missing parameter in model XML");
  }

  auto indep_source_num = settings->FirstChildElement("number-of-independent-sources");
  if (indep_source_num == NULL) {
    std::cerr << "The model xml misses the parameter number-of-independent-sources, "
            "child of the parameter model-settings." << std::endl;
    throw std::logic_error("Missing parameter in model XML");
  }

  return true;
}

//static bool IsValidAlgoXML(std::string path);
//static bool IsValidDataXML(std::string path);

bool InputsAssert::IsFileCorrect(std::string path, bool is_xml){
  return IsFilePathCorrect(path) && !IsFileEmpty(&path[0]) && (!is_xml || IsXMLValid(&path[0]));
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