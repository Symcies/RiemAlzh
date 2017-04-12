#include "ReadData_UT.h"

extern const std::string GV::TEST_DATA_DIR;

namespace test {
  void ReadData_UT::SetUp() {
    Test::SetUp();
  }

  TEST_F(ReadData_UT, read_obs) {
    std::string p = GV::TEST_DATA_DIR + "correct_real_multivariate_data_settings.xml";
    io::RealDataSettings data_settings(&p[0]);

    Observations obs = io::ReadData::ReadObservations(data_settings);
    obs.InitializeGlobalAttributes();

    ASSERT_FLOAT_EQ(obs.GetNumberOfSubjects(), 248);
    ASSERT_FLOAT_EQ(obs.GetTotalNumberOfObservations(), 1488);
    ASSERT_FLOAT_EQ(obs.GetTotalSumOfCognitiveScores(), 688.846);
    ASSERT_FLOAT_EQ(obs.GetTotalSumOfLandmarks(),0);

    //TODO: add part for each individual

  }

  TEST_F(ReadData_UT, open_kernel) {
    //To write when kernel data available
  }

}
