#ifndef MATRIX_TEST_MATRIX_SUBSCRIPT_H_
#define MATRIX_TEST_MATRIX_SUBSCRIPT_H_

#include <gtest/gtest.h>
#include "slab/matrix.h"

namespace slab {
TEST(MatrixTest, SubscriptByFortranStyle) {
  Matrix<int, 2> m{{01, 02, 03}, {11, 12, 13}};
  m(1, 2) = 99;

  EXPECT_EQ(1, m(0, 0));
  EXPECT_EQ(2, m(0, 1));
  EXPECT_EQ(3, m(0, 2));
  EXPECT_EQ(11, m(1, 0));
  EXPECT_EQ(12, m(1, 1));
  EXPECT_EQ(99, m(1, 2));
}

}  // namespace slab

#endif  // MATRIX_TEST_MATRIX_SUBSCRIPT_H_
