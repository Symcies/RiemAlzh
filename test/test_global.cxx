#include "test_global.h"

namespace GVT {
  extern const std::string BUILD_DIR = "/Users/clementine.fourrier/build";
  extern const std::string TEST_DIR = "/Users/clementine.fourrier/RiemAlzh/test/";

  extern const std::string TEST_DATA_DIR = TEST_DIR + "datasets/data/";
  extern const std::string TEST_MODEL_DIR = TEST_DIR + "datasets/models/";
  extern const std::string TEST_ALGO_DIR = TEST_DIR + "datasets/algorithm/";
  extern const std::string TEST_OUTPUTS_DIR = TEST_DIR + "datasets/outputs/";

  /// MODEL
  extern const std::string MULTIVAR_MODEL_CORRECT = TEST_MODEL_DIR + "multivariate_model_correct_settings.xml";
  extern const std::string MULTIVAR_MODEL_UNPARSABLE_XML = TEST_MODEL_DIR + "multivariate_model_wrong_settings_unparsable_xml.xml";
  extern const std::string MULTIVAR_MODEL_MISSING_PARAM = TEST_MODEL_DIR + "multivariate_model_wrong_settings_missing_param.xml";

  extern const std::string UNIVAR_MODEL_CORRECT = TEST_MODEL_DIR + "univariate_model_correct_settings.xml";
  extern const std::string UNIVAR_MODEL_UNPARSABLE_XML = TEST_MODEL_DIR + "univariate_model_wrong_settings_unparsable_xml.xml";
  extern const std::string UNIVAR_MODEL_MISSING_PARAM = TEST_MODEL_DIR + "univariate_model_wrong_settings_missing_param.xml";

  extern const std::string EMPTY_MODEL = TEST_MODEL_DIR + "empty_model_settings.xml" ;

  /// ALGORITHM
  extern const std::string ALGORITHM_CORRECT = TEST_ALGO_DIR + "algorithm_correct_settings.xml";
  extern const std::string ALGORITHM_INCORRECT_PARAM = TEST_ALGO_DIR + "algorithm_negative_iterations.xml";
  extern const std::string ALGORITHM_MISSING_PARAM = TEST_ALGO_DIR + "algorithm_missing_param.xml";
  extern const std::string ALGORITHM_UNPARSABLE_XML = TEST_ALGO_DIR + "algorithm_unparsable_xml.xml";
  extern const std::string EMPTY_ALGORITHM = TEST_ALGO_DIR + "empty_algorithm_settings.xml";

  /// DATA
  extern const std::string UNIVAR_DATA_CORRECT = TEST_DATA_DIR + "univariate_data_correct_real_settings.xml";
  extern const std::string MULTIVAR_DATA_CORRECT = TEST_DATA_DIR + "multivariate_data_correct_real_settings.xml";
  extern const std::string SIMULATED_DATA_CORRECT = TEST_DATA_DIR + "data_correct_simulated_settings.xml";

}
