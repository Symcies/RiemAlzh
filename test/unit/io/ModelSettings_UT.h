#pragma once

#include "gtest/gtest.h"
#include <iostream>
#include <fstream>
#include "ModelSettings.h"
#include "test/test_global.h"

namespace test {
  class ModelSettings_UT : public ::testing::Test {

    protected:
      virtual void SetUp();

  };
}
