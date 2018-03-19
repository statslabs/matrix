#include <gtest/gtest.h>
#include "test_matrix_construct.h"
#include "test_matrix_assign.h"
#include "test_matrix_subscript.h"
#include "test_matrix_opereration.h"

#include "test_blas.h"

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}