#pragma once

#include "gtest/gtest.h"
#include <fstream>
#include <cstdio>
#include "global.h"

namespace test {
  class ScalarModels_IT : public ::testing::Test {

    protected:
      virtual void SetUp();
  };
}