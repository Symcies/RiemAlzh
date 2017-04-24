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
//  PyRun_SimpleString("sys.path.append('/Users/clementine.fourrier/miniconda2/lib/python2.7/')");
  PyRun_SimpleString(("sys.path.append('"  + GV::PYTHON27_DIR + "lib-dynload/')").c_str());

  PySys_SetArgv(0, argv);

}

PythonUtils::~PythonUtils() {
  Py_Finalize();

}

void PythonUtils::CallPythonScript(std::string output_file_name) {
  PyObject *module_name, *module_obj, *module_contents_dict, *module_func, *module_arg, *output_path;

  // Nom du module
  std::string name = "PlotFinalOutput";

  // Convert the file name to a Python string.
  module_name = PyString_FromString(name.c_str());

  // Import the file as a Python module.
  module_obj = PyImport_Import(module_name);
  if (module_obj == nullptr) {
    PyErr_Print();
    std::cerr << "Fails to import the module.\n";
  }

  // Create a dictionary for the contents of the module.
  module_contents_dict = PyModule_GetDict(module_obj);

  // Get the plot function from the dictionary.
  module_func = PyDict_GetItemString(module_contents_dict, "Plot");

  // Create a Python tuple to hold the arguments to the method (None, in this case).
  module_arg = PyTuple_New(1);
  std::string output = GV::BUILD_DIR + output_file_name;
  output_path = PyString_FromString(output.c_str());
  PyTuple_SetItem(module_arg, 0, output_path);
  
  // Call the function with the arguments.
  PyObject_CallObject(module_func, module_arg);

}