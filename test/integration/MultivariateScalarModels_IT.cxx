#include "MultivariateScalarModels_IT.h"

#include <iostream>
#include <fstream>

#include "Algorithm.h"
#include "AlgorithmSettings.h"
#include "BlockedGibbsSampler.h"
#include "DataSettings.h"
#include "ModelSettings.h"
#include "MultivariateModel.h"
#include "UnivariateModel.h"
#include "Observations.h"


extern const std::string GV::BUILD_DIR;
extern const std::string GV::TEST_DIR;
extern const std::string GV::TEST_DATA_DIR;
extern const std::string GV::TEST_MODEL_DIR;
extern const std::string GV::TEST_ALGO_DIR;

namespace test {
  void MultivariateScalarModels_IT::SetUp() {
    Test::SetUp();
    std::remove(&(GV::TEST_DIR + "log_multivariate_file.txt")[0]);
  }

  void MultivariateScalarModels_IT::TearDown(){
    Test::TearDown();
    std::remove(&(GV::TEST_DIR + "log_multivariate_file.txt")[0]);
  }

 TEST_F(MultivariateScalarModels_IT, execution_of_multivariate_model_on_correct_real_dataset) {
    /// Load the file arguments
    io::ModelSettings     model_settings(&(GV::TEST_MODEL_DIR + "correct_multivariate_model_settings.xml")[0]);
    io::AlgorithmSettings algo_settings(&(GV::TEST_ALGO_DIR + "correct_algorithm_settings.xml")[0]);
    io::RealDataSettings  data_settings(&(GV::TEST_DATA_DIR + "correct_real_multivariate_data_settings.xml")[0]);
    /// Initialize the sampler
    std::shared_ptr<AbstractSampler> sampler = std::make_shared<BlockedGibbsSampler>();

    /// Initialize the model
    std::shared_ptr<AbstractModel> model;
    ASSERT_EQ(model_settings.GetType(), "Multivariate");
    model = std::make_shared<MultivariateModel>(model_settings);

    /// Initialize the data
    ASSERT_EQ(data_settings.IsReal(), true);
    Observations obs;
    obs = io::ReadData::ReadObservations(data_settings);
    obs.InitializeGlobalAttributes();
    ASSERT_FLOAT_EQ(obs.GetNumberOfSubjects(), 248);
    ASSERT_FLOAT_EQ(obs.GetTotalNumberOfObservations(), 1488);
    ASSERT_FLOAT_EQ(obs.GetTotalSumOfCognitiveScores(), 688.846);
    ASSERT_FLOAT_EQ(obs.GetTotalSumOfLandmarks(),0);

    /// Algorithm pipeline
    auto algo = std::make_shared<Algorithm>(algo_settings, model, sampler);
    algo->ComputeMCMCSAEM(obs);

    std::ifstream reference_file(GV::TEST_OUTPUTS_DIR + "ref_log_multivariate_real_file.txt");
    std::ifstream generated_file(GV::TEST_DIR + "log_multivariate_file.txt");
    float ref_value, gen_value;

    while(reference_file >> ref_value){
      generated_file >> gen_value;
      ASSERT_EQ(ref_value, gen_value);
    }
    std::remove(&(GV::TEST_DIR + "log_multivariate_file.txt")[0]);
  }

  TEST_F(MultivariateScalarModels_IT, execution_of_multivariate_model_on_correct_simulated_dataset) {
    /// Load the file arguments
    io::ModelSettings         model_settings(&(GV::TEST_MODEL_DIR + "correct_multivariate_model_settings.xml")[0]);
    io::AlgorithmSettings     algo_settings(&(GV::TEST_ALGO_DIR + "correct_algorithm_settings.xml")[0]);
    io::SimulatedDataSettings data_settings(&(GV::TEST_DATA_DIR + "correct_simulated_data_settings.xml")[0]);

    /// Initialize the sampler
    std::shared_ptr<AbstractSampler> sampler = std::make_shared<BlockedGibbsSampler>();

    /// Initialize the model
    std::shared_ptr<AbstractModel> model;
    ASSERT_EQ(model_settings.GetType(), "Multivariate");
    model = std::make_shared<MultivariateModel>(model_settings);

    /// Initialize the data
    ASSERT_EQ(data_settings.IsReal(), false);
    Observations obs;
    obs = model->SimulateData(data_settings, true);
    ASSERT_FLOAT_EQ(obs.GetNumberOfSubjects(), 150);

    /// Algorithm pipeline
    auto algo = std::make_shared<Algorithm>(algo_settings, model, sampler);
    algo->ComputeMCMCSAEM(obs);

    std::ifstream reference_file(GV::TEST_OUTPUTS_DIR + "ref_log_multivariate_simulated_file.txt");
    std::ifstream generated_file(GV::TEST_DIR + "log_multivariate_file.txt");
    float ref_value, gen_value;

    while(reference_file >> ref_value){
      generated_file >> gen_value;
      ASSERT_EQ(ref_value, gen_value);
    }
    std::remove(&(GV::TEST_DIR + "log_multivariate_file.txt")[0]);
  }

  TEST_F(MultivariateScalarModels_IT, multivariate_model_with_incorrect_params_model_settings) {
    /// Load the file arguments
    io::ModelSettings         model_settings(&(GV::TEST_MODEL_DIR + "incorrect_params_multivariate_model_settings.xml")[0]);
    io::AlgorithmSettings     algo_settings(&(GV::TEST_ALGO_DIR + "correct_algorithm_settings.xml")[0]);
    io::SimulatedDataSettings data_settings(&(GV::TEST_DATA_DIR + "correct_simulated_data_settings.xml")[0]);

    /// Initialize the sampler
    std::shared_ptr<AbstractSampler> sampler = std::make_shared<BlockedGibbsSampler>();

    /// Initialize the model
    std::shared_ptr<AbstractModel> model;
    ASSERT_EQ(model_settings.GetType(), "Multivariate");
    model = std::make_shared<MultivariateModel>(model_settings);

    /// Initialize the data
    ASSERT_EQ(data_settings.IsReal(), false);
    Observations obs;
    obs = model->SimulateData(data_settings, true);
    ASSERT_FLOAT_EQ(obs.GetNumberOfSubjects(), 150);

    /// Algorithm pipeline
    auto algo = std::make_shared<Algorithm>(algo_settings, model, sampler);
    algo->ComputeMCMCSAEM(obs);

    std::ifstream reference_file(GV::TEST_OUTPUTS_DIR + "ref_log_multivariate_simulated_file.txt");
    std::ifstream generated_file(GV::TEST_DIR + "log_multivariate_file.txt");
    float ref_value, gen_value;

    while(reference_file >> ref_value){
      generated_file >> gen_value;
      ASSERT_EQ(ref_value, gen_value);
    }
    std::remove(&(GV::TEST_DIR + "log_multivariate_file.txt")[0]);
  }
}
