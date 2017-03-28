#pragma once

#include "gtest/gtest.h"
#include <iostream>
#include <fstream>
#include "ReadData.h"

namespace test {
  class TestReadData : public ::testing::Test {

    protected:
      virtual void SetUp();

  };
}
