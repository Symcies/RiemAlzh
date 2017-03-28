#pragma once

#include "gtest/gtest.h"
#include <iostream>
#include <fstream>
#include "ModelSettings.h"

namespace test {
  class TestModelSettings : public ::testing::Test {

    protected:
      virtual void SetUp();

  };
}
