#include "TestModelSettings.h"

namespace test {
  void TestModelSettings::SetUp() {
    Test::SetUp();
  }

  TEST_F(TestModelSettings, construction) {
    const char * p1 = "/Users/clementine.fourrier/RiemAlzh/examples/scalar_models/MultivariateModel/model_settings.xml";
    io::ModelSettings     model_settings(p1);

    ASSERT_EQ(model_settings.GetType(), "Multivariate");
    ASSERT_EQ(model_settings.GetIndependentSourcesNumber(), 7);
    ASSERT_EQ(model_settings.GetInvertKernelPath(), "");
    ASSERT_EQ(model_settings.GetInterpolationKernelPath(), "");

  }

}
