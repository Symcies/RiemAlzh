#pragma once

#include "gtest/gtest.h"
#include <iostream>
#include <fstream>

namespace test {
  class TestScalarModel : public ::testing::Test {

    protected:
      virtual void SetUp();

  };
}