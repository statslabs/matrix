#include <gtest/gtest.h>
#include "test_matrix_ctor.h"
#include "test_matrix_assign.h"

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}