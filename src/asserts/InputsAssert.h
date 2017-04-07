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
  static bool IsFilePathCorrect(const char* path);
  static bool IsFilePathCorrect(std::string path);
  static bool IsXMLValid(const char* path);
  static std::string ToLowerCase(std::string str);
  static bool StringToBool(std::string str);

private:
  InputsAssert();
  InputsAssert(const InputsAssert&);
  InputsAssert& operator=(const InputsAssert&);

};