/* This class aims to provide all the function needed
 * to test the validity of the read inputs
 */

#pragma once

#include <fstream>
#include <iostream>
#include <string>

#include "tinyxml2.h"

#include "InputException.h"

class InputsAssert {
public:
  /// Basic file assertions
  static void IsFileCorrect(std::string path, bool is_xml);
  static void IsFilePathCorrect(std::string path);
  static void IsFileEmpty(std::string path);

  /// XML files assertions
  static void IsXMLValid(std::string path);
  static void IsValidModelXML(std::string path);
  static void IsValidAlgoXML(std::string path);
  static void IsValidSamplerXML(std::string path);
  static void IsValidDataXML(std::string path);
  static void IsValidSimulatedData(const tinyxml2::XMLElement * settings);
  static void IsValidRealData(const tinyxml2::XMLElement * settings);

  /// String modifications
  static std::string ToLowerCase(std::string str);
  static bool StringToBool(std::string str);


private:
  InputsAssert();
  InputsAssert(const InputsAssert&);
  InputsAssert& operator=(const InputsAssert&);

};