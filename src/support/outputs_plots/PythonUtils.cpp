#include "PythonUtils.h"
extern const std::string GV::SRC_DIR;
extern const std::string GV::PYTHON27_DIR;

PythonUtils::PythonUtils(char ** argv) {
  Py_SetPythonHome(strdup(GV::PYTHON27_DIR.c_str()));
  Py_SetProgramName("LongitudinaPython");

  Py_Initialize();
  PyRun_SimpleString("import sys");
  PyRun_SimpleString(("sys.path.append('" + GV::SRC_DIR + "support/outputs_plots/')").c_str());
  PyRun_SimpleString(("sys.path.append('" + GV::PYTHON27_DIR + "site-packages/matplotlib/')").c_str());
  PyRun_SimpleString(("sys.path.append('" + GV::PYTHON27_DIR + "site-packages/')").c_str());
  PyRun_SimpleString(("sys.path.append('"  + GV::PYTHON27_DIR + "lib-dynload/')").c_str());
  PyRun_SimpleString(("sys.path.append('"  + GV::PYTHON27_DIR + "lib-tk/')").c_str());

  PySys_SetArgv(0, argv);

  module_name_ = "PlotGraphs";

}

PythonUtils::~PythonUtils() {
  Py_Finalize();

}

void PythonUtils::PlotFinalOutput(std::string output_file_name, std::string type) {
  PyObject *module_name, *module_obj, *module_contents_dict, *module_func, *module_arg;

  // Convert the file name to a Python string.
  module_name = PyString_FromString(module_name_.c_str());

  // Import the file as a Python module.
  module_obj = PyImport_Import(module_name);
  if (module_obj == nullptr) {
    PyErr_Print();
    std::cerr << "Fails to import the module.\n";
  }

  // Create a dictionary for the contents of the module.
  module_contents_dict = PyModule_GetDict(module_obj);
  if (module_contents_dict == nullptr) {
    PyErr_Print();
    std::cerr << "Fails with module dictionnary.\n";
  }

  // Get the plot function from the dictionary.
  module_func = PyDict_GetItemString(module_contents_dict, "PlotFinalOutput");
  if (module_func == nullptr) {
    PyErr_Print();
    std::cerr << "Fails to import the function from the module.\n";
  }

  // Create a Python tuple to hold the arguments to the method (None, in this case).
  module_arg = PyTuple_New(2);
  std::string output = GV::BUILD_DIR + output_file_name;
  PyTuple_SetItem(module_arg, 0, PyString_FromString(output.c_str()));
  PyTuple_SetItem(module_arg, 1, PyString_FromString(type.c_str()));

  // Call the function with the arguments.
  PyObject_CallObject(module_func, module_arg);

}

void PythonUtils::PlotOutputWhileComputing(std::string output_file_name, int max_num_iter) {
  PyObject *module_name, *module_obj, *module_contents_dict, *module_func, *module_arg;

  // Convert the file name to a Python string.
  module_name = PyString_FromString(module_name_.c_str());

  // Import the file as a Python module.
  module_obj = PyImport_Import(module_name);
  if (module_obj == nullptr) {
    PyErr_Print();
    std::cerr << "Fails to import the module.\n";
  }

  // Create a dictionary for the contents of the module.
  module_contents_dict = PyModule_GetDict(module_obj);
  if (module_contents_dict == nullptr) {
    PyErr_Print();
    std::cerr << "Fails with module dictionnary.\n";
  }

  module_func = PyDict_GetItemString(module_contents_dict, "PlotOutputWhileComputing");
  if (module_func == nullptr) {
    PyErr_Print();
    std::cerr << "Fails to import the function from the module.\n";
  }

  // Create a Python tuple to hold the arguments to the method
  module_arg = PyTuple_New(1);
  std::string output = GV::BUILD_DIR + output_file_name;
  PyTuple_SetItem(module_arg, 0, PyString_FromString(output.c_str()));

  // Call the function with the arguments.
  PyObject_CallObject(module_func, module_arg);

}


void PythonUtils::PlotAllFinalOutputWithPatientData(std::string output_file_name, std::string type, Observations obs){
  PyObject *module_name, *module_obj, *module_contents_dict, *module_func, *module_arg;
  PyObject *obs_list, *observations, *key, *item;

  // Convert the file name to a Python string.
  module_name = PyString_FromString(module_name_.c_str());

  // Import the file as a Python module.
  module_obj = PyImport_Import(module_name);
  if (module_obj == nullptr) {
    PyErr_Print();
    std::cerr << "Fails to import the module.\n";
  }

  // Create a dictionary for the contents of the module.
  module_contents_dict = PyModule_GetDict(module_obj);
  if (module_contents_dict == nullptr) {
    PyErr_Print();
    std::cerr << "Fails with module dictionnary.\n";
  }

  // Get the plot function from the dictionary.
  module_func = PyDict_GetItemString(module_contents_dict, "PlotAndSelectPatientCurvesWithData");
  if (module_func == nullptr) {
    PyErr_Print();
    std::cerr << "Fails to import the function from the module.\n";
  }

  // Observations management
  obs_list = PyList_New(0);
  for(int i = 0; i< obs.GetNumberOfSubjects(); i++){
    observations = PyDict_New();
    for(int j = 0; j <obs.GetSubjectObservations(i).GetNumberOfTimePoints(); j++){
      key = PyFloat_FromDouble(obs.GetSubjectObservations(i).GetTimePoint(j));
      item = PyFloat_FromDouble(obs.GetSubjectObservations(i).GetCognitiveScore(j)[0]);
      PyDict_SetItem(observations, key, item);
    }
    PyList_Append(obs_list, observations);
  }

  // Create a Python tuple to hold the arguments to the method.
  module_arg = PyTuple_New(3);
  std::string output = GV::BUILD_DIR + output_file_name;
  PyTuple_SetItem(module_arg, 0, PyString_FromString(output.c_str()));
  PyTuple_SetItem(module_arg, 1, obs_list);
  PyTuple_SetItem(module_arg, 2, PyString_FromString(type.c_str()));

  // Call the function with the arguments.
  PyObject_CallObject(module_func, module_arg);

};