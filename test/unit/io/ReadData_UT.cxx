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

    ASSERT_FLOAT_EQ(obs.GetNumberOfSubjects(), 731);
    ASSERT_FLOAT_EQ(obs.GetTotalNumberOfObservations(), 3680);
    ASSERT_FLOAT_EQ(obs.GetTotalSumOfCognitiveScores(), 1048.1747); //
    ASSERT_FLOAT_EQ(obs.GetTotalSumOfLandmarks(),0);

    //TODO: add part for each individual

  }


}
