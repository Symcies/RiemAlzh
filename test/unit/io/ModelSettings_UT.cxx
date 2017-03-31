#include "ModelSettings_UT.h"

namespace test {
  void ModelSettings_UT::SetUp() {
    Test::SetUp();
  }

  TEST_F(ModelSettings_UT, construction_multivariate_model) {
    const char * p1 = "/Users/clementine.fourrier/RiemAlzh/test/datasets/models/correct_multivariate_model_settings.xml";
    io::ModelSettings     model_settings(p1);

    ASSERT_EQ(model_settings.GetType(), "Multivariate");
    ASSERT_EQ(model_settings.GetIndependentSourcesNumber(), 7);
    ASSERT_EQ(model_settings.GetInvertKernelPath(), "");
    ASSERT_EQ(model_settings.GetInterpolationKernelPath(), "");

  }

}
