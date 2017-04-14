#pragma once

#include "gtest/gtest.h"

#include <cstdio>
#include <iostream>
#include <fstream>

#include "global.h"
#include "test/test_global.h"

#include "fit.h"

#include "Algorithm.h"
#include "AlgorithmSettings.h"
#include "BlockedGibbsSampler.h"
#include "DataSettings.h"
#include "ModelSettings.h"
#include "MultivariateModel.h"
#include "UnivariateModel.h"
#include "Observations.h"

namespace test {
  class MultivariateScalarModels_IT : public ::testing::Test {

    protected:
      virtual void SetUp();
      virtual void TearDown();
  };
}
