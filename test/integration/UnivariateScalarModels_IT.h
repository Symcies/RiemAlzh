#pragma once

#include "gtest/gtest.h"
#include <fstream>
#include <cstdio>
#include "test/test_global.h"

namespace test {
  class UnivariateScalarModels_IT : public ::testing::Test {

    protected:
      virtual void SetUp();
      virtual void TearDown();
  };
}
