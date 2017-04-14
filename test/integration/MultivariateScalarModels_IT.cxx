#include "MultivariateScalarModels_IT.h"

extern const std::string GVT::BUILD_DIR;
extern const std::string GVT::TEST_DIR;
extern const std::string GVT::TEST_DATA_DIR;
extern const std::string GVT::TEST_MODEL_DIR;
extern const std::string GVT::TEST_ALGO_DIR;
extern const std::string GVT::MULTIVAR_MODEL_CORRECT;
extern const std::string GVT::MULTIVAR_MODEL_UNPARSABLE_XML;
extern const std::string GVT::MULTIVAR_MODEL_MISSING_PARAM;
extern const std::string GVT::UNIVAR_MODEL_CORRECT;
extern const std::string GVT::UNIVAR_MODEL_UNPARSABLE_XML;
extern const std::string GVT::UNIVAR_MODEL_MISSING_PARAM;
extern const std::string GVT::EMPTY_MODEL;

namespace test {
  void MultivariateScalarModels_IT::SetUp() {
    Test::SetUp();
    std::remove((GVT::TEST_DIR + "log_multivariate_file.txt").c_str());
  }

  void MultivariateScalarModels_IT::TearDown(){
    Test::TearDown();
    std::remove((GVT::TEST_DIR + "log_multivariate_file.txt").c_str());
  }

  /// CORRECT EXECUTION TESTS

  TEST_F(MultivariateScalarModels_IT, correct_real_dataset) {/*
    try {
  //TODO: will have to change in order not to use strdup
      char* params[] = {"Longitudina", "fit", strdup(GVT::MULTIVAR_MODEL_CORRECT.c_str()), strdup(GVT::ALGORITHM_CORRECT.c_str()),
                        strdup(GVT::MULTIVAR_DATA_CORRECT.c_str()),strdup(GVT::SAMPLER_CORRECT.c_str())};
      fit(6, params);
    } catch(InputException exception) {
      FAIL() << "Received exception " << exception.what();
    } catch(std::exception exception) {
      FAIL() << "Received exception " << exception.what();
    }

    /// End result assertions
    std::ifstream reference_file(GVT::TEST_OUTPUTS_DIR + "ref_log_multivariate_real_file.txt");
    std::ifstream generated_file(GVT::TEST_DIR + "log_multivariate_file.txt");
    float ref_value, gen_value;

    while(reference_file >> ref_value){
      generated_file >> gen_value;
      ASSERT_EQ(ref_value, gen_value);
    }
    std::remove((GVT::TEST_DIR + "log_multivariate_file.txt").c_str());*/
  }

  TEST_F(MultivariateScalarModels_IT, correct_simulated_dataset) {
    //TODO: replace with the pipeline when existing
    try {
      //TODO: will have to change in order not to use strdup
      /// Load the file arguments
      io::ModelSettings model_settings((GVT::MULTIVAR_MODEL_CORRECT).c_str());
      io::AlgorithmSettings algo_settings((GVT::ALGORITHM_CORRECT).c_str());
      io::SimulatedDataSettings data_settings((GVT::SIMULATED_DATA_CORRECT).c_str());
      io::SamplerSettings sampler_settings((GVT::SAMPLER_CORRECT).c_str());

      /// Initialize the model
      std::shared_ptr <AbstractModel> model;
      model = std::make_shared<MultivariateModel>(model_settings);
std::cout << "2" << std::endl;

      /// Initialize the data
      Observations obs;
      obs = model->SimulateData(data_settings, true);
std::cout << "3" << std::endl;

      /// Algorithm pipeline
      auto algo = std::make_shared<Algorithm>(algo_settings);
std::cout << "4" << std::endl;
      algo->SetModel(model);
std::cout << "5" << std::endl;
      algo->AddSamplers(sampler_settings);
std::cout << "6" << std::endl;
      algo->ComputeMCMCSAEM(obs);
std::cout << "7" << std::endl;
    }  catch(InputException exception) {
      FAIL() << "Received exception " << exception.what();
    } catch(std::exception exception) {
      FAIL() << "Received exception " << exception.what();
    }

    std::ifstream reference_file((GVT::TEST_OUTPUTS_DIR + "ref_log_multivariate_simulated_file.txt").c_str());
    std::ifstream generated_file((GVT::TEST_DIR + "log_multivariate_file.txt").c_str());
    float ref_value, gen_value;

    /// Final result assertions
    while(reference_file >> ref_value){
      generated_file >> gen_value;
      ASSERT_EQ(ref_value, gen_value);
    }
    std::remove((GVT::TEST_DIR + "log_multivariate_file.txt").c_str());
  }

/* This area might be redundant with the unit failre tests
  /// EXECUTION PROBLEMS TEST IN MODEL SETTINGS

  TEST_F(MultivariateScalarModels_IT, missing_params_model_settings) {
    /// Load the file arguments
    bool error_detected = false;
    try {
      //TODO: replace with the pipeline when existing
      io::ModelSettings model_settings((GVT::MULTIVAR_MODEL_MISSING_PARAM).c_str());
      io::AlgorithmSettings algo_settings((GVT::ALGORITHM_CORRECT).c_str());
      io::RealDataSettings data_settings((GVT::MULTIVAR_DATA_CORRECT).c_str());
    } catch(InputException exception){
      ASSERT_STREQ(exception.what(), "The model xml misses the parameter type, child of the parameter model-settings.");
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }

  TEST_F(MultivariateScalarModels_IT, incorrect_xml_model_settings) {
  /// Load the file arguments
    bool error_detected = false;
    try {
      //TODO: replace with the pipeline when existing
      io::ModelSettings model_settings((GVT::MULTIVAR_MODEL_UNPARSABLE_XML).c_str());
      io::AlgorithmSettings algo_settings((GVT::ALGORITHM_CORRECT).c_str());
      io::RealDataSettings data_settings((GVT::MULTIVAR_DATA_CORRECT).c_str());
    } catch(InputException exception){
      std::string error_message = "The XML of file " + GVT::MULTIVAR_MODEL_UNPARSABLE_XML + " is not correct.";
      ASSERT_STREQ(exception.what(), error_message.c_str());
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }


  TEST_F(MultivariateScalarModels_IT, empty_xml_model_settings) {
  /// Load the file arguments
    bool error_detected = false;
    try {
    //TODO: replace with the pipeline when existing
      io::ModelSettings model_settings((GVT::EMPTY_MODEL).c_str());
      io::AlgorithmSettings algo_settings((GVT::ALGORITHM_CORRECT).c_str());
      io::RealDataSettings data_settings((GVT::MULTIVAR_DATA_CORRECT).c_str());
    } catch(InputException exception){
      std::string error_message =  "The file at '" + GVT::EMPTY_MODEL + "' is empty.";
      ASSERT_STREQ(exception.what(), error_message.c_str());
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }


  TEST_F(MultivariateScalarModels_IT, incorrect_path_model_settings) {
  /// Load the file arguments
  bool error_detected = false;
  try {
    //TODO: replace with the pipeline when existing
      io::ModelSettings model_settings("this/is/not/a/path");
      io::AlgorithmSettings algo_settings((GVT::ALGORITHM_CORRECT).c_str());
      io::RealDataSettings data_settings((GVT::MULTIVAR_DATA_CORRECT).c_str());
    } catch(InputException exception){
      ASSERT_STREQ(exception.what(), "The file path 'this/is/not/a/path' is incorrect.");
      error_detected = true;
    } catch(std::exception exception){
      FAIL() << "Exception thrown was not of type InputException. Was " << exception.what();
    }
    ASSERT_EQ(error_detected, true);
  }

*/
}
