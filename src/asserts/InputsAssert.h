/* This class aims to privde all the function needed to test the validity of the read inputs
 *
 */

#pragma once

#include <fstream>
#include <iostream>
#include <string>

#include "tinyxml2.h"

class InputsAssert {
public:
  /// Basic file assertions
  static bool IsFileCorrect(std::string path, bool is_xml);
  static bool IsFilePathCorrect(std::string path);
  static bool IsFileEmpty(std::string path);

  /// XML files assertions
  static bool IsXMLValid(std::string path);
  static bool IsValidModelXML(std::string path);
  static bool IsValidAlgoXML(std::string path);
  static bool IsValidDataXML(std::string path);

  /// String modifications
  static std::string ToLowerCase(std::string str);
  static bool StringToBool(std::string str);


private:
  InputsAssert();
  InputsAssert(const InputsAssert&);
  InputsAssert& operator=(const InputsAssert&);

};