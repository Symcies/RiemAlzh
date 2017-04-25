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

  void PlotFinalOutput(std::string type);
  void PlotOutputWhileComputing(std::string output_file_name, std::string type);

};
