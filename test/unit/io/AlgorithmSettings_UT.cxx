#include "AlgorithmSettings_UT.h"

extern const std::string GV::TEST_ALGO_DIR;

namespace test {
  void AlgorithmSettings_UT::SetUp() {
    Test::SetUp();
  }

  TEST_F(AlgorithmSettings_UT, construction_algorithm) {
    std::string path = GV::TEST_ALGO_DIR + "correct_algorithm_settings.xml";
    io::AlgorithmSettings algo_settings(&path[0]);

    ASSERT_EQ(algo_settings.GetMaximumNumberOfIterations(), 1000);
    ASSERT_EQ(algo_settings.GetNumberOfBurnInIterations(),  25000);
    ASSERT_EQ(algo_settings.GetOutputDisplayIteration(),    0);
    ASSERT_EQ(algo_settings.GetDataSaveIteration(),         100);

  }

}
