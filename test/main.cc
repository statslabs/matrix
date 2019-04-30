#include <gtest/gtest.h>
#include "test_blas.h"
#include "test_construction_and_assignment.h"
#include "test_matrix_functions.h"
#include "test_matrix_opereration.h"
#include "test_matrix_subscript.h"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
