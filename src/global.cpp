//Global variables management
#include "global.h"

namespace GV{ //for global variables
  extern const std::string BUILD_DIR = "/Users/clementine.fourrier/build";
  extern const std::string TEST_DIR = "/Users/clementine.fourrier/RiemAlzh/test/";

  //extern const std::string BUILD_DIR = "/Users/igor.koval/Documents/Work";
  //extern const std::string TEST_DIR = "/Users/igor.koval/Documents/Work/RiemAlzh/test/";
  extern const std::string TEST_DATA_DIR = TEST_DIR + "datasets/data/";
  extern const std::string TEST_MODEL_DIR = TEST_DIR + "datasets/models/";
  extern const std::string TEST_ALGO_DIR = TEST_DIR + "datasets/algorithm/";
  extern const std::string TEST_OUTPUTS_DIR = TEST_DIR + "datasets/outputs/";
  extern bool        TEST_RUN = true;
}

