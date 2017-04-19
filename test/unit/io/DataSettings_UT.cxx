#include "DataSettings_UT.h"

extern const std::string GVT::TEST_DATA_DIR;
extern const std::string GVT::TEST_MODEL_DIR;
extern const std::string GVT::TEST_ALGO_DIR;

namespace test {
  void DataSettings_UT::SetUp() {
    Test::SetUp();
  }

  TEST_F(DataSettings_UT, construction_from_real_data) {
    std::string p = GVT::MULTIVAR_DATA_CORRECT;
    io::RealDataSettings data_settings(&p[0]);

    ASSERT_EQ(
      data_settings.GetPathToGroup(),
      GVT::TEST_DATA_DIR + "data_files_multivariate/group.csv"
    );
    ASSERT_EQ(
      data_settings.GetPathToTimepoints(),
      GVT::TEST_DATA_DIR + "data_files_multivariate/X.csv"
    );
    ASSERT_EQ(
      data_settings.GetPathToCognitiveScores(),
      GVT::TEST_DATA_DIR + "data_files_multivariate/Y.csv"
    );
    ASSERT_EQ(data_settings.GetCognitiveScoresDimension(), 4);
    ASSERT_EQ(data_settings.GetPathToLandmarks(), "");

    ASSERT_EQ(data_settings.GetLandmarksDimension(), 0);

    ASSERT_EQ(data_settings.IsReal(), true);
    ASSERT_EQ(data_settings.LandmarkPresence(), false);
    ASSERT_EQ(data_settings.CognitiveScoresPresence(), true);

  }

  TEST_F(DataSettings_UT, construction_from_simulated_data) {
    std::string p = GVT::SIMULATED_DATA_CORRECT;
    io::SimulatedDataSettings data_settings(&p[0]);
    
    ASSERT_EQ(data_settings.IsReal(), false);
    ASSERT_EQ(data_settings.GetNumberOfSimulatedSubjects(), 150);
    ASSERT_EQ(data_settings.GetMinimumNumberOfObservations(), 4);
    ASSERT_EQ(data_settings.GetMaximumNumberOfObservations(), 6);
  }

/// EXECUTION PROBLEMS TEST IN MODEL SETTINGS

  TEST_F(DataSettings_UT, missing_params_data_settings) {
    /// Load the file arguments
    bool error_detected = false;
    try {
      io::RealDataSettings data_settings((GVT::MULTIVAR_DATA_MISSING_PARAM).c_str());
    } catch(InputException exception){
      ASSERT_STREQ(exception.what(), "The data xml misses the parameter group-file, child of the parameter real-data.");
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }

  TEST_F(DataSettings_UT, incorrect_xml_data_settings) {
    /// Load the file arguments
    bool error_detected = false;
    try {
      io::RealDataSettings data_settings((GVT::MULTIVAR_DATA_UNPARSABLE_XML).c_str());
    } catch(InputException exception){
      std::string error_message = "The XML of file " + GVT::MULTIVAR_DATA_UNPARSABLE_XML + " is not correct.";
      ASSERT_STREQ(exception.what(), error_message.c_str());
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }


  TEST_F(DataSettings_UT, empty_xml_data_settings) {
    /// Load the file arguments
    bool error_detected = false;
    try {
      io::RealDataSettings data_settings((GVT::EMPTY_MODEL).c_str());
    } catch(InputException exception){
      std::string error_message =  "The file at '" + GVT::EMPTY_MODEL + "' is empty.";
      ASSERT_STREQ(exception.what(), error_message.c_str());
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }


  TEST_F(DataSettings_UT, incorrect_path_data_settings) {
    /// Load the file arguments
    bool error_detected = false;
    try {
      io::RealDataSettings data_settings("this/is/not/a/path");
    } catch(InputException exception){
      ASSERT_STREQ(exception.what(), "The file path 'this/is/not/a/path' is incorrect.");
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }

  TEST_F(DataSettings_UT, incorrect_type_in_data_settings) {
    /// Load the file arguments
    bool error_detected = false;
    try {
      io::RealDataSettings data_settings((GVT::MULTIVAR_DATA_INCORRECT_PARAM_TYPE).c_str());
    } catch(InputException exception){
      std::string error_message = "test should be 'true' or 'false'.";
      ASSERT_STREQ(exception.what(), error_message.c_str());
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }

}
