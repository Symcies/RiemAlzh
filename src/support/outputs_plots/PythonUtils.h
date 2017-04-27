#pragma once

#include <iostream>
#include <sstream>
#include "stdio.h"
#include <string>
#include <Python.h>
#include "global.h"
#include "IndividualObservations.h"
#include "Observations.h"

class PythonUtils {
public:
  PythonUtils(char ** argv);
  ~PythonUtils();

  void PlotFinalOutput(std::string output_file_name, std::string type);
  void PlotOneFinalOutputWithPatientData(std::string output_file_name, std::string type, IndividualObservations obs, int obs_num);
  void PlotAllFinalOutputWithPatientData(std::string output_file_name, std::string type, Observations obs);
  void PlotOutputWhileComputing(std::string output_file_name, int max_num_iter);

private:
  std::string module_name_;

};
