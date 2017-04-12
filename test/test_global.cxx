#include "test_global.h"

namespace GVT {
  extern const std::string BUILD_DIR = "/Users/clementine.fourrier/build";
  extern const std::string TEST_DIR = "/Users/clementine.fourrier/RiemAlzh/test/";

  extern const std::string TEST_DATA_DIR = TEST_DIR + "datasets/data/";
  extern const std::string TEST_MODEL_DIR = TEST_DIR + "datasets/models/";
  extern const std::string TEST_ALGO_DIR = TEST_DIR + "datasets/algorithm/";
  extern const std::string TEST_OUTPUTS_DIR = TEST_DIR + "datasets/outputs/";

  extern const std::string MULTIVAR_MODEL_CORRECT = TEST_MODEL_DIR + "multivariate_model_correct_settings.xml";
  extern const std::string MULTIVAR_MODEL_UNPARSABLE_XML = TEST_MODEL_DIR + "multivariate_model_wrong_settings_unparsable_xml.xml";
  extern const std::string MULTIVAR_MODEL_MISSING_PARAM = TEST_MODEL_DIR + "multivariate_model_wrong_settings_missing_param.xml";

  extern const std::string UNIVAR_MODEL_CORRECT = TEST_MODEL_DIR + "univariate_model_correct_settings.xml";
  extern const std::string UNIVAR_MODEL_UNPARSABLE_XML = TEST_MODEL_DIR + "univariate_model_wrong_settings_unparsable_xml.xml";
  extern const std::string UNIVAR_MODEL_MISSING_PARAM = TEST_MODEL_DIR + "univariate_model_wrong_settings_missing_param.xml";

  extern const std::string EMPTY_MODEL = TEST_MODEL_DIR + "empty_model_settings.xml" ;
}
