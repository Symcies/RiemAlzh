#include "gtest/gtest.h"
#include "global.h"
#include "test_global.h"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  GV::TEST_RUN = true;

  return RUN_ALL_TESTS();

}
