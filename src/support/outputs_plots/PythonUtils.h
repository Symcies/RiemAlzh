#pragma once

#include <iostream>
#include <sstream>
#include "stdio.h"
#include <string>
#include <Python.h>
#include "global.h"

class PythonUtils {
public:
  PythonUtils(char ** argv);
  ~PythonUtils();

  void CallPythonScript(std::string output_file_name);

};
