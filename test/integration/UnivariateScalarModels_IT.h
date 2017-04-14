#pragma once

#include "gtest/gtest.h"

#include <cstdio>
#include <fstream>
#include <iostream>

#include "fit.h"

#include "test/test_global.h"

#include "Algorithm.h"
#include "AlgorithmSettings.h"
#include "BlockedGibbsSampler.h"
#include "DataSettings.h"
#include "ModelSettings.h"
#include "MultivariateModel.h"
#include "UnivariateModel.h"
#include "Observations.h"

namespace test {
  class UnivariateScalarModels_IT : public ::testing::Test {

    protected:
      virtual void SetUp();
      virtual void TearDown();
  };
}
