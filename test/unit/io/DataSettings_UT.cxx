#include "DataSettings_UT.h"

extern const std::string GV::TEST_DATA_DIR;
extern const std::string GV::TEST_MODEL_DIR;
extern const std::string GV::TEST_ALGO_DIR;

namespace test {
  void DataSettings_UT::SetUp() {
    Test::SetUp();
  }

  TEST_F(DataSettings_UT, construction_from_real_data) {
    std::string p = GV::TEST_DATA_DIR + "correct_real_multivariate_data_settings.xml";
    io::RealDataSettings data_settings(&p[0]);

    ASSERT_EQ(
      data_settings.GetPathToGroup(),
      GV::TEST_DATA_DIR + "data_files_multivariate/group.csv"
    );
    ASSERT_EQ(
      data_settings.GetPathToTimepoints(),
      GV::TEST_DATA_DIR + "data_files_multivariate/X.csv"
    );
    ASSERT_EQ(
      data_settings.GetPathToCognitiveScores(),
      GV::TEST_DATA_DIR + "data_files_multivariate/Y.csv"
    );
    ASSERT_EQ(
      data_settings.GetCognitiveScoresDimension(),
      4
    );
    ASSERT_EQ(
      data_settings.GetPathToLandmarks(),
      ""
    );

    ASSERT_EQ(
      data_settings.GetLandmarksDimension(),
      0
    );

    ASSERT_EQ(
      data_settings.IsReal(),
      true
    );
    ASSERT_EQ(
      data_settings.LandmarkPresence(),
      false
    );
    ASSERT_EQ(
      data_settings.CognitiveScoresPresence(),
      true
    );

  }

  TEST_F(DataSettings_UT, construction_from_simulated_data) {
    std::string p = GV::TEST_DATA_DIR + "correct_simulated_data_settings.xml";
    io::SimulatedDataSettings data_settings(&p[0]);


    ASSERT_EQ(
      data_settings.IsReal(),
      false
    );
    ASSERT_EQ(
      data_settings.GetNumberOfSimulatedSubjects(),
      150
    );
    ASSERT_EQ(
      data_settings.GetMinimumNumberOfObservations(),
      4
    );
    ASSERT_EQ(
      data_settings.GetMaximumNumberOfObservations(),
      6
    );
  }

}
