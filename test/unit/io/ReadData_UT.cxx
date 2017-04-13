#include "ReadData_UT.h"

extern const std::string GVT::TEST_DATA_DIR;

namespace test {
  void ReadData_UT::SetUp() {
    Test::SetUp();
  }

  TEST_F(ReadData_UT, read_obs) {
    std::string p = GVT::MULTIVAR_DATA_CORRECT;
    io::RealDataSettings data_settings(p.c_str());

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
