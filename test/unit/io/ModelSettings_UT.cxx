#include "ModelSettings_UT.h"

extern const std::string GVT::TEST_MODEL_DIR;

namespace test {
  void ModelSettings_UT::SetUp() {
    Test::SetUp();
  }

/// WORKING TEST
  TEST_F(ModelSettings_UT, construction_multivariate_model) {
    io::ModelSettings model_settings(GVT::MULTIVAR_MODEL_CORRECT.c_str());

    ASSERT_EQ(model_settings.GetType(), "Multivariate");
    ASSERT_EQ(model_settings.GetIndependentSourcesNumber(), 7);
    ASSERT_EQ(model_settings.GetInvertKernelPath(), "");
    ASSERT_EQ(model_settings.GetInterpolationKernelPath(), "");
  }

/// EXECUTION PROBLEMS TEST IN MODEL SETTINGS

  TEST_F(ModelSettings_UT, missing_params_model_settings) {
    /// Load the file arguments
    bool error_detected = false;
    try {
      io::ModelSettings model_settings((GVT::MULTIVAR_MODEL_MISSING_PARAM).c_str());
    } catch(InputException exception){
      ASSERT_STREQ(exception.what(), "The model xml misses the parameter type, child of the parameter model-settings.");
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }

  TEST_F(ModelSettings_UT, missing_variables_model_settings) {
    /// Load the file arguments
    bool error_detected = false;
    try {
      io::ModelSettings model_settings((GVT::MULTIVAR_MODEL_MISSING_VAR).c_str());
    } catch(InputException exception){
      ASSERT_STREQ(exception.what(), "The model xml misses the parameter noise");
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }

  TEST_F(ModelSettings_UT, incorrect_xml_model_settings) {
    /// Load the file arguments
    bool error_detected = false;
    try {
      io::ModelSettings model_settings((GVT::MULTIVAR_MODEL_UNPARSABLE_XML).c_str());
    } catch(InputException exception){
      std::string error_message = "The XML of file " + GVT::MULTIVAR_MODEL_UNPARSABLE_XML + " is not correct.";
      ASSERT_STREQ(exception.what(), error_message.c_str());
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }


  TEST_F(ModelSettings_UT, empty_xml_model_settings) {
    /// Load the file arguments
    bool error_detected = false;
    try {
      io::ModelSettings model_settings((GVT::EMPTY_MODEL).c_str());
    } catch(InputException exception){
      std::string error_message =  "The file at '" + GVT::EMPTY_MODEL + "' is empty.";
      ASSERT_STREQ(exception.what(), error_message.c_str());
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }


  TEST_F(ModelSettings_UT, incorrect_path_model_settings) {
    /// Load the file arguments
    bool error_detected = false;
    try {
      io::ModelSettings model_settings("this/is/not/a/path");
    } catch(InputException exception){
      ASSERT_STREQ(exception.what(), "The file path 'this/is/not/a/path' is incorrect.");
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }

  TEST_F(ModelSettings_UT, incorrect_type_in_model_settings) {
    /// Load the file arguments
    bool error_detected = false;
    try {
      io::ModelSettings model_settings((GVT::MULTIVAR_MODEL_INCORRECT_PARAM_TYPE).c_str());
    } catch(InputException exception){
      std::string error_message = "The parameter mean is not a double.";
      ASSERT_STREQ(exception.what(), error_message.c_str());
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }


}
