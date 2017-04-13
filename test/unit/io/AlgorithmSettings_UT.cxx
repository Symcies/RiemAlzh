#include "AlgorithmSettings_UT.h"

extern const std::string GVT::TEST_ALGO_DIR;

namespace test {
  void AlgorithmSettings_UT::SetUp() {
    Test::SetUp();
  }

  TEST_F(AlgorithmSettings_UT, construction_algorithm) {
    std::string path = GVT::ALGORITHM_CORRECT;
    io::AlgorithmSettings algo_settings(path.c_str());

    ASSERT_EQ(algo_settings.GetMaximumNumberOfIterations(), 1000);
    ASSERT_EQ(algo_settings.GetNumberOfBurnInIterations(),  25000);
    ASSERT_EQ(algo_settings.GetOutputDisplayIteration(),    0);
    ASSERT_EQ(algo_settings.GetDataSaveIteration(),         100);

  }


/// EXECUTION PROBLEMS TEST IN ALGORITHM SETTINGS

  TEST_F(AlgorithmSettings_UT, missing_params_algo_settings) {
    /// Load the file arguments
    bool error_detected = false;
    try {
      io::AlgorithmSettings algo_settings((GVT::ALGORITHM_MISSING_PARAM).c_str());
    } catch(InputException exception){
      ASSERT_STREQ(exception.what(), "The algorithm xml misses the parameter max-iterations, child of the parameter algorithm-settings.");
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }

  TEST_F(AlgorithmSettings_UT, incorrect_xml_algo_settings) {
  /// Load the file arguments
    bool error_detected = false;
    try {
      io::AlgorithmSettings algo_settings((GVT::ALGORITHM_UNPARSABLE_XML).c_str());
    } catch(InputException exception){
      std::string error_message = "The XML of file " + GVT::ALGORITHM_UNPARSABLE_XML + " is not correct.";
      ASSERT_STREQ(exception.what(), error_message.c_str());
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }


  TEST_F(AlgorithmSettings_UT, negative_unsigned_int_algo_settings) {
    /// Load the file arguments
    bool error_detected = false;
    try {
      io::AlgorithmSettings algo_settings((GVT::ALGORITHM_INCORRECT_PARAM).c_str());
    } catch(InputException exception){
      ASSERT_STREQ(exception.what(), "All iterations values must be superior or equal to 0.");
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }

  TEST_F(AlgorithmSettings_UT, empty_xml_algo_settings) {
    /// Load the file arguments
    bool error_detected = false;
    try {
      io::AlgorithmSettings algo_settings((GVT::EMPTY_ALGORITHM).c_str());
    } catch(InputException exception){
      std::string error_message =  "The file at '" + GVT::EMPTY_ALGORITHM + "' is empty.";
      ASSERT_STREQ(exception.what(), error_message.c_str());
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }


  TEST_F(AlgorithmSettings_UT, incorrect_path_algo_settings) {
    /// Load the file arguments
    bool error_detected = false;
    try {
      io::AlgorithmSettings algo_settings("this/is/not/a/path");
    } catch(InputException exception){
      ASSERT_STREQ(exception.what(), "The file path 'this/is/not/a/path' is incorrect.");
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }


}
