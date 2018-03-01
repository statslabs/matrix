#ifndef MATRIX_TEST_MATRIX_ASSIGN_H_
#define MATRIX_TEST_MATRIX_ASSIGN_H_

#include <gtest/gtest.h>
#include "slab/matrix.h"

namespace slab {
TEST(MatrixAssignTest, AssignFromMatrixRef) {
  Matrix<double, 1> m1{1, 2, 3};
  auto mr1 = m1(slice(1));                    // MatrixRef<double, 1>
  Matrix<double, 1> m1sub = mr1;

  EXPECT_EQ(2, m1sub(0));
  EXPECT_EQ(3, m1sub(1));

  Matrix<double, 2> m2{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  auto mr2 = m2(slice(1), slice(1));          // MatrixRef<double, 2>
  Matrix<double, 2> m2sub = mr2;

  EXPECT_EQ(5, m2sub(0, 0));
  EXPECT_EQ(6, m2sub(0, 1));
  EXPECT_EQ(8, m2sub(1, 0));
  EXPECT_EQ(9, m2sub(1, 1));

  Matrix<double, 3>
      m3{{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
         {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
         {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}};
  auto mr3 = m3(slice(1), slice(1), slice(1)); // MatrixRef<double, 3>
  Matrix<double, 3> m3sub = mr3;

  EXPECT_EQ(5, m3sub(0, 0, 0));
  EXPECT_EQ(6, m3sub(0, 0, 1));
  EXPECT_EQ(8, m3sub(0, 1, 0));
  EXPECT_EQ(9, m3sub(0, 1, 1));
  EXPECT_EQ(5, m3sub(1, 0, 0));
  EXPECT_EQ(6, m3sub(1, 0, 1));
  EXPECT_EQ(8, m3sub(1, 1, 0));
  EXPECT_EQ(9, m3sub(1, 1, 1));
}

TEST(MatrixAssignTest, AssignFromMatrixInitializer) {
  Matrix<double, 0> m0 = {1};

  EXPECT_EQ(m0(), 1);

  Matrix<double, 1> m1 = {1, 2};

  EXPECT_EQ(m1(0), 1);
  EXPECT_EQ(m1(1), 2);

  Matrix<double, 2> m2 = {{1, 2}, {3, 4}};

  EXPECT_EQ(m2(0, 0), 1);
  EXPECT_EQ(m2(0, 1), 2);
  EXPECT_EQ(m2(1, 0), 3);
  EXPECT_EQ(m2(1, 1), 4);

  Matrix<double, 3> m3 = {{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}};

  EXPECT_EQ(m3(0, 0, 0), 1);
  EXPECT_EQ(m3(0, 0, 1), 2);
  EXPECT_EQ(m3(0, 1, 0), 3);
  EXPECT_EQ(m3(0, 1, 1), 4);
  EXPECT_EQ(m3(1, 0, 0), 5);
  EXPECT_EQ(m3(1, 0, 1), 6);
  EXPECT_EQ(m3(1, 1, 0), 7);
  EXPECT_EQ(m3(1, 1, 1), 8);
}


}

#endif //MATRIX_TEST_MATRIX_ASSIGN_H
