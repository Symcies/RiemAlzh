#pragma once

#include "gtest/gtest.h"
#include <iostream>
#include <fstream>
#include "ReadData.h"
#include "RealDataSettings.h"
#include "test/test_global.h"

namespace test {
  class ReadData_UT : public ::testing::Test {

    protected:
      virtual void SetUp();

  };
}
