//
// Created by Yi Pan (Institute of Cancer and Genomic Sciences) on 12/03/2018.
//

#ifndef MATRIX_TEST_MATRIX_OPERERATION_H
#define MATRIX_TEST_MATRIX_OPERERATION_H

#include <gtest/gtest.h>
#include "slab/matrix.h"

namespace slab {

TEST(MatrixOperationTest, Diag_Vec) {
  vec v = {1, 2, 3};
  mat m = diag(v);

  EXPECT_EQ(1, m(0, 0));
  EXPECT_EQ(0, m(0, 1));
  EXPECT_EQ(0, m(0, 2));
  EXPECT_EQ(0, m(1, 0));
  EXPECT_EQ(2, m(1, 1));
  EXPECT_EQ(0, m(1, 2));
  EXPECT_EQ(0, m(2, 0));
  EXPECT_EQ(0, m(2, 1));
  EXPECT_EQ(3, m(2, 2));
}

TEST(MatrixOperationTest, Transpose_Mat) {
  mat m1 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  mat m2 = transpose(m1);

  EXPECT_EQ(1, m2(0, 0));
  EXPECT_EQ(4, m2(0, 1));
  EXPECT_EQ(7, m2(0, 2));
  EXPECT_EQ(2, m2(1, 0));
  EXPECT_EQ(5, m2(1, 1));
  EXPECT_EQ(8, m2(1, 2));
  EXPECT_EQ(3, m2(2, 0));
  EXPECT_EQ(6, m2(2, 1));
  EXPECT_EQ(9, m2(2, 2));
}

TEST(MatrixOperationTest, Matmul_Mat_Vec) {
  mat m1 = {{8, 4, 7}, {3, 5, 1}, {1, 3, 2}};
  vec v1 = {-1, 2, 1};
  vec v2 = {3, 2, 3};

  vec res = matmul(m1, v1) + v2;

  EXPECT_EQ(10, res(0));
  EXPECT_EQ(10, res(1));
  EXPECT_EQ(10, res(2));
}

TEST(MatrixOperationTest, Matmul_Mat_Mat) {
  mat m1 = {{1, 2, 3}, {4, 5, 6}};
  mat m2 = {{7, 8}, {9, 10}, {11, 12}};
  mat res = matmul(m1, m2);

  EXPECT_EQ(58, res(0, 0));
  EXPECT_EQ(64, res(0, 1));
  EXPECT_EQ(139, res(1, 0));
  EXPECT_EQ(154, res(1, 1));
}

TEST(MatrixOperationTest, Matmul_N) {
  mat m1 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  mat m2 = m1;
  mat m3 = m1;

  mat res = matmul_n(m1, m2, m3);

  EXPECT_EQ(468, res(0, 0));
  EXPECT_EQ(576, res(0, 1));
  EXPECT_EQ(684, res(0, 2));

  EXPECT_EQ(1062, res(1, 0));
  EXPECT_EQ(1305, res(1, 1));
  EXPECT_EQ(1548, res(1, 2));

  EXPECT_EQ(1656, res(2, 0));
  EXPECT_EQ(2034, res(2, 1));
  EXPECT_EQ(2412, res(2, 2));
}

TEST(MatrixOperationTest, Sum) {
  mat m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

  double res = sum(m);

  EXPECT_EQ(45, res);
}

TEST(MatrixOperationTest, Exp) {
  mat m = zeros<mat>(3, 3);
  mat res = exp(m);

  EXPECT_EQ(1, res(0, 0));
  EXPECT_EQ(1, res(0, 1));
  EXPECT_EQ(1, res(0, 2));
  EXPECT_EQ(1, res(1, 0));
  EXPECT_EQ(1, res(1, 1));
  EXPECT_EQ(1, res(1, 2));
  EXPECT_EQ(1, res(2, 0));
  EXPECT_EQ(1, res(2, 1));
  EXPECT_EQ(1, res(2, 2));
}

TEST(MatrixOperationTest, Log) {
  mat m = ones<mat>(3, 3);
  mat res = log(m);

  EXPECT_EQ(0, res(0, 0));
  EXPECT_EQ(0, res(0, 1));
  EXPECT_EQ(0, res(0, 2));
  EXPECT_EQ(0, res(1, 0));
  EXPECT_EQ(0, res(1, 1));
  EXPECT_EQ(0, res(1, 2));
  EXPECT_EQ(0, res(2, 0));
  EXPECT_EQ(0, res(2, 1));
  EXPECT_EQ(0, res(2, 2));
}

TEST(MatrixOperationTest, Pow) {
  mat m = ones<mat>(3, 3) * 2.0;
  mat res = pow(m, -1);

  EXPECT_EQ(0.5, res(0, 0));
  EXPECT_EQ(0.5, res(0, 1));
  EXPECT_EQ(0.5, res(0, 2));
  EXPECT_EQ(0.5, res(1, 0));
  EXPECT_EQ(0.5, res(1, 1));
  EXPECT_EQ(0.5, res(1, 2));
  EXPECT_EQ(0.5, res(2, 0));
  EXPECT_EQ(0.5, res(2, 1));
  EXPECT_EQ(0.5, res(2, 2));
}

}  // namespace slab

#endif  // MATRIX_TEST_MATRIX_OPERERATION_H
