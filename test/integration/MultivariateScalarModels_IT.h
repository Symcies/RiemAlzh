#pragma once

#include "gtest/gtest.h"
#include <fstream>
#include <cstdio>
#include "global.h"
#include "test/test_global.h"

namespace test {
  class MultivariateScalarModels_IT : public ::testing::Test {

    protected:
      virtual void SetUp();
      virtual void TearDown();
  };
}
