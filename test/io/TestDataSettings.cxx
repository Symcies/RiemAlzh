#include "TestDataSettings.h"

namespace test {
  void TestDataSettings::SetUp() {
    Test::SetUp();
  }

  TEST_F(TestDataSettings, construction_from_real_data) {
    const char * p = "/Users/clementine.fourrier/RiemAlzh/examples/scalar_models/MultivariateModel/data_settings.xml";
    io::DataSettings data_settings(p);

    ASSERT_EQ(
      data_settings.GetPathToGroup(),
      "/Users/clementine.fourrier/RiemAlzh/examples/scalar_models/MultivariateModel/data/group.csv"
    );
    ASSERT_EQ(
      data_settings.GetPathToTimepoints(),
      "/Users/clementine.fourrier/RiemAlzh/examples/scalar_models/MultivariateModel/data/X.csv"
    );
    ASSERT_EQ(
      data_settings.GetPathToCognitiveScores(),
      "/Users/clementine.fourrier/RiemAlzh/examples/scalar_models/MultivariateModel/data/Y.csv"
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
      -1
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


    ASSERT_EQ(
      data_settings.GetNumberOfSimulatedSubjects(),
      -1
    );
    ASSERT_EQ(
      data_settings.GetMinimumNumberOfObservations(),
      -1
    );
    ASSERT_EQ(
      data_settings.GetMaximumNumberOfObservations(),
      -1
    );
  }

  TEST_F(TestDataSettings, construction_from_false_data) {
    const char * p = "/Users/clementine.fourrier/RiemAlzh/examples/scalar_models/MultivariateModel/data_settings0.xml";
    io::DataSettings data_settings(p);

    ASSERT_EQ(
      data_settings.GetPathToGroup(),
      ""
    );
    ASSERT_EQ(
      data_settings.GetPathToTimepoints(),
      ""
    );
    ASSERT_EQ(
      data_settings.GetPathToCognitiveScores(),
      ""
    );
    ASSERT_EQ(
      data_settings.GetCognitiveScoresDimension(),
      -1
    );
    ASSERT_EQ(
      data_settings.GetPathToLandmarks(),
      ""
    );

    ASSERT_EQ(
      data_settings.GetLandmarksDimension(),
      -1
    );

    ASSERT_EQ(
      data_settings.IsReal(),
      false
    );
    ASSERT_EQ(
      data_settings.LandmarkPresence(),
      false
    );
    ASSERT_EQ(
      data_settings.CognitiveScoresPresence(),
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
