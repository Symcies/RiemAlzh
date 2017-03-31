#include "ScalarModels_IT.h"

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

namespace test {
  void ScalarModels_IT::SetUp() {
    Test::SetUp();
  }

  // TEST_F(ScalarModels_IT, execution_of_multivariate_model_on_correct_real_dataset) {
    /// Load the file arguments
  //   io::ModelSettings     model_settings(&(TEST_DIR + MODEL_DIR + "correct_multivariate_model_settings.xml")[0]);
  //   io::AlgorithmSettings algo_settings(&(TEST_DIR + ALGO_DIR + "correct_algorithm_settings.xml")[0]);
  //   io::DataSettings      data_settings(&(TEST_DIR + DATA_DIR + "correct_real_data_settings.xml")[0]);
  //
  //   /// Initialize the sampler
  //   std::shared_ptr<AbstractSampler> sampler = std::make_shared<BlockedGibbsSampler>();
  //
  //   /// Initialize the model
  //   std::shared_ptr<AbstractModel> model;
  //   ASSERT_EQ(model_settings.GetType(), "Multivariate");
  //   model = std::make_shared<MultivariateModel>(model_settings);
  //
  //   /// Initialize the data
  //   ASSERT_EQ(data_settings.IsReal(), true);
  //   Observations obs;
  //   obs = io::ReadData::ReadObservations(data_settings);
  //   obs.InitializeGlobalAttributes();
  //   ASSERT_FLOAT_EQ(obs.GetNumberOfSubjects(), 248);
  //   ASSERT_FLOAT_EQ(obs.GetTotalNumberOfObservations(), 1488);
  //   ASSERT_FLOAT_EQ(obs.GetTotalSumOfCognitiveScores(), 688.846);
  //   ASSERT_FLOAT_EQ(obs.GetTotalSumOfLandmarks(),0);
  //
  //   /// Algorithm pipeline
  //   auto algo = std::make_shared<Algorithm>(algo_settings);
  //   algo->SetModel(model);
  //   algo->SetSampler(sampler);
  //   algo->ComputeMCMCSAEM(obs);
  //
  //   std::string result[] =
  //     {"Iteration n: 0 - noise: 0.34233 - G: 0.11986 - T0: 75.1052 - Var(Tau): 0.183834 - Ksi: 0.130696 - Var(Ksi): 0.000465735",
  //     "Iteration n: 100 - noise: 0.181954 - G: 0.120081 - T0: 77.906 - Var(Tau): 12.2058 - Ksi: 0.12972 - Var(Ksi): 0.000650701",
  //     "Iteration n: 200 - noise: 0.147939 - G: 0.120133 - T0: 78.7398 - Var(Tau): 16.216 - Ksi: 0.131095 - Var(Ksi): 0.000715965",
  //     "Iteration n: 300 - noise: 0.113185 - G: 0.120042 - T0: 79.4669 - Var(Tau): 21.6711 - Ksi: 0.127058 - Var(Ksi): 0.000884629",
  //     "Iteration n: 400 - noise: 0.0990816 - G: 0.119984 - T0: 79.9915 - Var(Tau): 23.9557 - Ksi: 0.123335 - Var(Ksi): 0.00070051",
  //     "Iteration n: 500 - noise: 0.0833372 - G: 0.120179 - T0: 79.9951 - Var(Tau): 26.3214 - Ksi: 0.122554 - Var(Ksi): 0.000595611",
  //     "Iteration n: 600 - noise: 0.0839312 - G: 0.120466 - T0: 80.4003 - Var(Tau): 26.2095 - Ksi: 0.115875 - Var(Ksi): 0.000619247",
  //     "Iteration n: 700 - noise: 0.0769561 - G: 0.120237 - T0: 80.4113 - Var(Tau): 28.3012 - Ksi: 0.112834 - Var(Ksi): 0.000600446",
  //     "Iteration n: 800 - noise: 0.0741597 - G: 0.120237 - T0: 80.2624 - Var(Tau): 27.9319 - Ksi: 0.107817 - Var(Ksi): 0.000509121",
  //     "Iteration n: 900 - noise: 0.0742886 - G: 0.120127 - T0: 80.2902 - Var(Tau): 30.7549 - Ksi: 0.101566 - Var(Ksi): 0.0006735"};
  //
  //   std::string line;
  //   std::ifstream is("log_file.txt");
  //   int iter = 0;
  //   while(std::getline(is, line), iter++){
  //     ASSERT_EQ(line, result[iter]);
  //   }
  // }

  TEST_F(ScalarModels_IT, execution_of_univariate_model_on_correct_real_dataset) {
    /// Load the file arguments
    io::ModelSettings     model_settings(&(TEST_DIR + MODEL_DIR + "correct_univariate_model_settings.xml")[0]);
    io::AlgorithmSettings algo_settings(&(TEST_DIR + ALGO_DIR + "correct_algorithm_settings.xml")[0]);
    io::DataSettings      data_settings(&(TEST_DIR + DATA_DIR + "correct_real_data_settings.xml")[0]);

    /// Initialize the sampler
    std::shared_ptr<AbstractSampler> sampler = std::make_shared<BlockedGibbsSampler>();

    /// Initialize the model
    std::shared_ptr<AbstractModel> model;
    ASSERT_EQ(model_settings.GetType(), "Univariate");
    model = std::make_shared<UnivariateModel>(model_settings);

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
    auto algo = std::make_shared<Algorithm>(algo_settings);
    algo->SetModel(model);
    algo->SetSampler(sampler);
    algo->ComputeMCMCSAEM(obs);
    //
    // std::string result[] =
    //   {"Iteration n: 0 - noise: 0.34233 - G: 0.11986 - T0: 75.1052 - Var(Tau): 0.183834 - Ksi: 0.130696 - Var(Ksi): 0.000465735",
    //   "Iteration n: 100 - noise: 0.181954 - G: 0.120081 - T0: 77.906 - Var(Tau): 12.2058 - Ksi: 0.12972 - Var(Ksi): 0.000650701",
    //   "Iteration n: 200 - noise: 0.147939 - G: 0.120133 - T0: 78.7398 - Var(Tau): 16.216 - Ksi: 0.131095 - Var(Ksi): 0.000715965",
    //   "Iteration n: 300 - noise: 0.113185 - G: 0.120042 - T0: 79.4669 - Var(Tau): 21.6711 - Ksi: 0.127058 - Var(Ksi): 0.000884629",
    //   "Iteration n: 400 - noise: 0.0990816 - G: 0.119984 - T0: 79.9915 - Var(Tau): 23.9557 - Ksi: 0.123335 - Var(Ksi): 0.00070051",
    //   "Iteration n: 500 - noise: 0.0833372 - G: 0.120179 - T0: 79.9951 - Var(Tau): 26.3214 - Ksi: 0.122554 - Var(Ksi): 0.000595611",
    //   "Iteration n: 600 - noise: 0.0839312 - G: 0.120466 - T0: 80.4003 - Var(Tau): 26.2095 - Ksi: 0.115875 - Var(Ksi): 0.000619247",
    //   "Iteration n: 700 - noise: 0.0769561 - G: 0.120237 - T0: 80.4113 - Var(Tau): 28.3012 - Ksi: 0.112834 - Var(Ksi): 0.000600446",
    //   "Iteration n: 800 - noise: 0.0741597 - G: 0.120237 - T0: 80.2624 - Var(Tau): 27.9319 - Ksi: 0.107817 - Var(Ksi): 0.000509121",
    //   "Iteration n: 900 - noise: 0.0742886 - G: 0.120127 - T0: 80.2902 - Var(Tau): 30.7549 - Ksi: 0.101566 - Var(Ksi): 0.0006735"};
    //
    // std::string line;
    // std::ifstream is("log_file.txt");
    // int iter = 0;
    // while(std::getline(is, line), iter++){
    //   ASSERT_EQ(line, result[iter]);
    // }
  }

  TEST_F(ScalarModels_IT, execution_of_multivariate_model_on_correct_simulated_dataset) {
    /// Load the file arguments
    io::ModelSettings     model_settings(&(TEST_DIR + MODEL_DIR + "correct_multivariate_model_settings.xml")[0]);
    io::AlgorithmSettings algo_settings(&(TEST_DIR + ALGO_DIR + "correct_algorithm_settings.xml")[0]);
    io::DataSettings      data_settings(&(TEST_DIR + DATA_DIR + "correct_simulated_data_settings.xml")[0]);

    /// Initialize the sampler
    std::shared_ptr<AbstractSampler> sampler = std::make_shared<BlockedGibbsSampler>();

    /// Initialize the model
    std::shared_ptr<AbstractModel> model;
    ASSERT_EQ(model_settings.GetType(), "Multivariate");
    model = std::make_shared<MultivariateModel>(model_settings);

    /// Initialize the data
    ASSERT_EQ(data_settings.IsReal(), false);
    Observations obs;
    obs = io::ReadData::ReadObservations(data_settings);
    obs.InitializeGlobalAttributes();
    // ASSERT_FLOAT_EQ(obs.GetNumberOfSubjects(), 150);
    // ASSERT_FLOAT_EQ(obs.GetTotalNumberOfObservations(), 1488);
    // ASSERT_FLOAT_EQ(obs.GetTotalSumOfCognitiveScores(), 688.846);
    // ASSERT_FLOAT_EQ(obs.GetTotalSumOfLandmarks(),0);

    /// Algorithm pipeline
    auto algo = std::make_shared<Algorithm>(algo_settings);
    algo->SetModel(model);
    algo->SetSampler(sampler);
    algo->ComputeMCMCSAEM(obs);

    std::string result[] =
      {"Iteration n: 0 - noise: 0.34233 - G: 0.11986 - T0: 75.1052 - Var(Tau): 0.183834 - Ksi: 0.130696 - Var(Ksi): 0.000465735",
      "Iteration n: 100 - noise: 0.181954 - G: 0.120081 - T0: 77.906 - Var(Tau): 12.2058 - Ksi: 0.12972 - Var(Ksi): 0.000650701",
      "Iteration n: 200 - noise: 0.147939 - G: 0.120133 - T0: 78.7398 - Var(Tau): 16.216 - Ksi: 0.131095 - Var(Ksi): 0.000715965",
      "Iteration n: 300 - noise: 0.113185 - G: 0.120042 - T0: 79.4669 - Var(Tau): 21.6711 - Ksi: 0.127058 - Var(Ksi): 0.000884629",
      "Iteration n: 400 - noise: 0.0990816 - G: 0.119984 - T0: 79.9915 - Var(Tau): 23.9557 - Ksi: 0.123335 - Var(Ksi): 0.00070051",
      "Iteration n: 500 - noise: 0.0833372 - G: 0.120179 - T0: 79.9951 - Var(Tau): 26.3214 - Ksi: 0.122554 - Var(Ksi): 0.000595611",
      "Iteration n: 600 - noise: 0.0839312 - G: 0.120466 - T0: 80.4003 - Var(Tau): 26.2095 - Ksi: 0.115875 - Var(Ksi): 0.000619247",
      "Iteration n: 700 - noise: 0.0769561 - G: 0.120237 - T0: 80.4113 - Var(Tau): 28.3012 - Ksi: 0.112834 - Var(Ksi): 0.000600446",
      "Iteration n: 800 - noise: 0.0741597 - G: 0.120237 - T0: 80.2624 - Var(Tau): 27.9319 - Ksi: 0.107817 - Var(Ksi): 0.000509121",
      "Iteration n: 900 - noise: 0.0742886 - G: 0.120127 - T0: 80.2902 - Var(Tau): 30.7549 - Ksi: 0.101566 - Var(Ksi): 0.0006735"};

    std::string line;
    std::ifstream is("log_file.txt");
    int iter = 0;
    while(std::getline(is, line), iter++){
      ASSERT_EQ(line, result[iter]);
    }
  }

}
