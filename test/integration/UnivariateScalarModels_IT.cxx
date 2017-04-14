#include "UnivariateScalarModels_IT.h"

extern const std::string GVT::BUILD_DIR;
extern const std::string GVT::TEST_DIR;
extern const std::string GVT::TEST_DATA_DIR;
extern const std::string GVT::TEST_MODEL_DIR;
extern const std::string GVT::TEST_ALGO_DIR;

namespace test {
  void UnivariateScalarModels_IT::SetUp() {
    Test::SetUp();
    std::remove((GVT::TEST_DIR + "log_univariate_file.txt").c_str());
  }

  void UnivariateScalarModels_IT::TearDown(){
    Test::TearDown();
    std::remove((GVT::TEST_DIR + "log_univariate_file.txt").c_str());
  }

  TEST_F(UnivariateScalarModels_IT, execution_of_univariate_model_on_correct_real_dataset) {
    /// Load the file arguments
    char* params[] = {"Longitudina", "fit", strdup(GVT::UNIVAR_MODEL_CORRECT.c_str()), strdup(GVT::ALGORITHM_CORRECT.c_str()),
                    strdup(GVT::UNIVAR_DATA_CORRECT.c_str()),strdup(GVT::SAMPLER_CORRECT.c_str())};
    fit(6, params);

    std::ifstream reference_file(GVT::TEST_OUTPUTS_DIR + "ref_log_univariate_real_file.txt");
    std::ifstream generated_file(GVT::TEST_DIR + "log_univariate_file.txt");
    float ref_value, gen_value;

    while(reference_file >> ref_value){
      generated_file >> gen_value;
      ASSERT_EQ(ref_value, gen_value);
    }
    std::remove((GVT::TEST_DIR + "log_univariate_file.txt").c_str());

  }

  TEST_F(UnivariateScalarModels_IT, execution_of_univariate_model_on_correct_simulated_dataset) {
    char* params[] = {"Longitudina", "fit", strdup(GVT::UNIVAR_MODEL_CORRECT.c_str()), strdup(GVT::ALGORITHM_CORRECT.c_str()),
                  strdup(GVT::SIMULATED_DATA_CORRECT.c_str()),strdup(GVT::SAMPLER_CORRECT.c_str())};
    fit(6, params);

    std::ifstream reference_file(GVT::TEST_OUTPUTS_DIR + "ref_log_univariate_simulated_file.txt");
    std::ifstream generated_file(GVT::TEST_DIR + "log_univariate_file.txt");
    float ref_value, gen_value;

    while(reference_file >> ref_value){
      generated_file >> gen_value;
      ASSERT_EQ(ref_value, gen_value);
    }
    std::remove((GVT::TEST_DIR + "log_univariate_file.txt").c_str());
  }

}
