#include "TestScalarModel.h"


namespace test {
  void TestScalarModel::SetUp() {
    Test::SetUp();
  }

  TEST_F(TestScalarModel, Test_google_test) {
    ASSERT_EQ(2+2, 4);
  }
}
