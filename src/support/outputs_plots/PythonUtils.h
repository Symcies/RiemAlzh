#pragma once

#include <iostream>
#include <sstream>
#include "stdio.h"
#include <string>
#include <Python.h>

class PythonUtils {
public:
  PythonUtils(char ** argv);
  ~PythonUtils();

  void CallPythonScript();

};
