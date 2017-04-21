#include "PythonUtils.h"

PythonUtils::PythonUtils(char ** argv) {
  Py_SetPythonHome("/Users/clementine.fourrier/miniconda2/lib/python2.7/");
  Py_SetProgramName("LongitudinaPython");
  Py_Initialize();
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("sys.path.append('/Users/clementine.fourrier/RiemAlzh/src/support/outputs_plots/')");
  PyRun_SimpleString("sys.path.append('/Users/clementine.fourrier/miniconda2/lib/python2.7/site-packages/numpy/')");
  PySys_SetArgv(0, argv);

}

PythonUtils::~PythonUtils() {
  Py_Finalize();

}

void PythonUtils::CallPythonScript() {
  std::string module_name = "PlotFinalOutput";

  //FILE * file;
  //file = fopen(file_path.c_str(), std::string("r").c_str());

  std::cout << "Before run?" << std::endl;
  //PyRun_AnyFile(file, file_name.c_str());

  PyObject *pName, *pModule, *pDict, *pFunc, *pArgs, *pValue;

// Convert the file name to a Python string.
  pName = PyString_FromString(module_name.c_str());

// Import the file as a Python module.
  pModule = PyImport_ImportModule("numpy");
  //pModule = PyImport_Import(pName);
  if (pModule == nullptr) {
    PyErr_Print();
    std::cerr << "Fails to import the module.\n";
  }
  std::cout << "name of module : " << pModule->ob_type->tp_name << std::endl;

// Create a dictionary for the contents of the module.
  pDict = PyModule_GetDict(pModule);

// Get the add method from the dictionary.
  pFunc = PyDict_GetItemString(pDict, "Plot");

// Create a Python tuple to hold the arguments to the method.
  pArgs = PyTuple_New(2);

  /*
// Convert 2 to a Python integer.
  pValue = PyInt_FromLong(2);

// Set the Python int as the first and second arguments to the method.
  PyTuple_SetItem(pArgs, 0, pValue);
  PyTuple_SetItem(pArgs, 1, pValue);*/

// Call the function with the arguments.

  PyObject* pResult = PyObject_CallObject(pFunc, pArgs);
  std::cout << "After run?" << std::endl;

}