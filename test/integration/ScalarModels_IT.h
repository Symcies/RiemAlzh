#pragma once

#include "gtest/gtest.h"
#include <fstream>

namespace test {
  class ScalarModels_IT : public ::testing::Test {

    protected:
      virtual void SetUp();
      std::string TEST_DIR = "/Users/clementine.fourrier/RiemAlzh/test/";
      std::string DATA_DIR = "datasets/data/";
      std::string MODEL_DIR = "datasets/models/";
      std::string ALGO_DIR = "datasets/algorithm/";

  };
}
