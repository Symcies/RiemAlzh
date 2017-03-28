#include "TestAlgorithmSettings.h"

namespace test {
  void TestAlgorithmSettings::SetUp() {
    Test::SetUp();
  }

  TEST_F(TestAlgorithmSettings, construction) {
    std::string path = "/Users/clementine.fourrier/RiemAlzh/examples/scalar_models/MultivariateModel/algorithm_settings.xml";
    io::AlgorithmSettings algo_settings(&path[0]);

    ASSERT_EQ(algo_settings.GetMaximumNumberOfIterations(), 1000);
    ASSERT_EQ(algo_settings.GetNumberOfBurnInIterations(),  25000);
    ASSERT_EQ(algo_settings.GetOutputDisplayIteration(),    0);
    ASSERT_EQ(algo_settings.GetDataSaveIteration(),         100);

  }

}
