#pragma once

#include "gtest/gtest.h"
#include <iostream>
#include <fstream>

#include "Algorithm.h"
#include "BlockedGibbsSampler.h"
#include "MultivariateModel.h"

namespace test {
  class TestAlgorithm : public ::testing::Test {

    protected:
      virtual void SetUp();
      std::shared_ptr<Algorithm> algo_;

  };
}
