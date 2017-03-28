#include "TestReadData.h"

namespace test {
  void TestReadData::SetUp() {
    Test::SetUp();
  }

  TEST_F(TestReadData, read_obs) {
    const char * p = "/Users/clementine.fourrier/RiemAlzh/examples/scalar_models/MultivariateModel/data_settings.xml";
    io::DataSettings data_settings(p);

    Observations obs = io::ReadData::ReadObservations(data_settings);
    obs.InitializeGlobalAttributes();

    ASSERT_FLOAT_EQ(obs.GetNumberOfSubjects(), 248);
    ASSERT_FLOAT_EQ(obs.GetTotalNumberOfObservations(), 1488);
    ASSERT_FLOAT_EQ(obs.GetTotalSumOfCognitiveScores(), 688.846);
    ASSERT_FLOAT_EQ(obs.GetTotalSumOfLandmarks(),0);

    //TODO: add part for each individual

  }

  TEST_F(TestReadData, open_kernel) {
    //To write when kernel data available
  }

}
