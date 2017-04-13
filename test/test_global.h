#pragma once
#include <string>

namespace GVT{ //for global variables test
  extern const std::string BUILD_DIR;
  extern const std::string TEST_DIR;

  extern const std::string TEST_DATA_DIR;
  extern const std::string TEST_MODEL_DIR;
  extern const std::string TEST_ALGO_DIR;
  extern const std::string TEST_OUTPUTS_DIR;

  extern const std::string MULTIVAR_MODEL_CORRECT;
  extern const std::string MULTIVAR_MODEL_UNPARSABLE_XML;
  extern const std::string MULTIVAR_MODEL_MISSING_PARAM;

  extern const std::string UNIVAR_MODEL_CORRECT;
  extern const std::string UNIVAR_MODEL_UNPARSABLE_XML;
  extern const std::string UNIVAR_MODEL_MISSING_PARAM;

  extern const std::string EMPTY_MODEL;


  extern const std::string ALGORITHM_CORRECT;
  extern const std::string ALGORITHM_INCORRECT_PARAM;
  extern const std::string ALGORITHM_MISSING_PARAM ;
  extern const std::string ALGORITHM_UNPARSABLE_XML;
  extern const std::string EMPTY_ALGORITHM;


  extern const std::string UNIVAR_DATA_CORRECT;
  extern const std::string MULTIVAR_DATA_CORRECT;
  extern const std::string SIMULATED_DATA_CORRECT;

  extern const std::string SAMPLER_CORRECT;

}
