#pragma once

#include "gtest/gtest.h"
#include <iostream>
#include <fstream>
#include "DataSettings.h"

namespace test {
  class TestDataSettings : public ::testing::Test {

    protected:
      virtual void SetUp();

  };
}
