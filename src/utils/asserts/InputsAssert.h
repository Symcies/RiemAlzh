/* This class aims to provide all the function needed
 * to test the validity of the read inputs
 */

#pragma once

#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <unordered_set>

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
  static void AllVariablesPresentInModelXML(std::string type, std::unordered_set<std::string> names);
  static void IsValidAlgoXML(std::string path);
  static void IsValidSamplerXML(std::string path);
  static void IsValidDataXML(std::string path);
  static void IsValidSimulatedData(const tinyxml2::XMLElement * settings);
  static void IsValidRealData(const tinyxml2::XMLElement * settings);

  /// File path regex checks
  static void IsCsvFilePath(std::string path);
  static void IsXMLFilePath(std::string path);

  /// String modifications
  static std::string ToLowerCase(std::string str);
  static bool StringToBool(std::string str);

  /// Input type checks
  static void CheckElemContentIsInt(const tinyxml2::XMLElement* elem);
  static void CheckElemContentIsString(const tinyxml2::XMLElement* elem);
  static void CheckElemContentIsBool(const tinyxml2::XMLElement* elem);
  static void CheckElemContentIsDouble(const tinyxml2::XMLElement* elem);

private:
  InputsAssert();
  InputsAssert(const InputsAssert&);
  InputsAssert& operator=(const InputsAssert&);

};