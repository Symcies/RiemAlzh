#pragma once

#include "gtest/gtest.h"
#include <iostream>
#include <fstream>
#include "ReadData.h"
#include "global.h"

namespace test {
  class ReadData_UT : public ::testing::Test {

    protected:
      virtual void SetUp();

  };
}
