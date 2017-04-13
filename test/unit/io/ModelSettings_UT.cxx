#include "ModelSettings_UT.h"

extern const std::string GV::TEST_MODEL_DIR;

namespace test {
  void ModelSettings_UT::SetUp() {
    Test::SetUp();
  }

  TEST_F(ModelSettings_UT, construction_multivariate_model) {
    std::string p1 = GV::TEST_MODEL_DIR + "correct_multivariate_model_settings.xml";
    io::ModelSettings     model_settings(&p1[0]);

    ASSERT_EQ(model_settings.GetType(), "Multivariate");
    ASSERT_EQ(model_settings.GetIndependentSourcesNumber(), 7);
    ASSERT_EQ(model_settings.GetInvertKernelPath(), "");
    ASSERT_EQ(model_settings.GetInterpolationKernelPath(), "");

  }

}
