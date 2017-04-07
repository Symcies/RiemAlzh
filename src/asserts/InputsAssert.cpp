#include "InputsAssert.h"


bool InputsAssert::IsFilePathCorrect(const char * path){
  std::ifstream test(path);
  if (test) {
    return true;
  }
  std::cerr << "The file path " << path << " is incorrect." << std::endl;
  return false;
}

bool InputsAssert::IsFilePathCorrect(std::string path){
  std::ifstream test(path);
  if (test) {
    return true;
  }
  std::cerr << "The file path " << path << " is incorrect." << std::endl;
  return false;
}

bool InputsAssert::IsXMLValid(const char * path){
  tinyxml2::XMLDocument file;
  file.LoadFile(path);

  if (file.Error()){
    std::cerr << "The XML of file " << path << " is not correct." <<std::endl;
  }
  return !file.Error();
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