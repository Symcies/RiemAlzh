#pragma once

#include "gtest/gtest.h"
#include <iostream>
#include <fstream>
#include "DataSettings.h"
#include "global.h"

namespace test {
  class DataSettings_UT : public ::testing::Test {

    protected:
      virtual void SetUp();

  };
}
