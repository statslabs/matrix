//
// Created by Yi Pan (Institute of Cancer and Genomic Sciences) on 12/03/2018.
//

#ifndef MATRIX_TEST_MATRIX_OPERERATION_H
#define MATRIX_TEST_MATRIX_OPERERATION_H

#include <gtest/gtest.h>
#include "slab/matrix.h"

namespace slab {

TEST(MatrixOperationTest, MatTranspose) {
  mat m1 = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9}
  };
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

TEST(MatrixOperationTest, IntMatVecProd) {
  imat m1 = {
      {8, 4, 7},
      {3, 5, 1},
      {1, 3, 2}
  };
  ivec v1 = {-1, 2, 1};
  ivec v2 = {3, 2, 3};

  ivec res = matmul(m1, v1) + v2;

  EXPECT_EQ(10, res(0));
  EXPECT_EQ(10, res(1));
  EXPECT_EQ(10, res(2));
}

TEST(MatrixOperationTest, FloatMatVecProd) {
  fmat m1 = {
      {8, 4, 7},
      {3, 5, 1},
      {1, 3, 2}
  };
  fvec v1 = {-1, 2, 1};
  fvec v2 = {3, 2, 3};

  fvec res = matmul(m1, v1) + v2;

  EXPECT_EQ(10, res(0));
  EXPECT_EQ(10, res(1));
  EXPECT_EQ(10, res(2));
}

TEST(MatrixOperationTest, DoubleMatVecProd) {
  mat m1 = {
      {8, 4, 7},
      {3, 5, 1},
      {1, 3, 2}
  };
  vec v1 = {-1, 2, 1};
  vec v2 = {3, 2, 3};

  vec res = matmul(m1, v1) + v2;

  EXPECT_EQ(10, res(0));
  EXPECT_EQ(10, res(1));
  EXPECT_EQ(10, res(2));
}

TEST(MatrixOperationTest, IntMatMatProd) {
  imat m1 = {
      {1, 2, 3},
      {4, 5, 6}
  };
  imat m2 = {
      {7, 8},
      {9, 10},
      {11, 12}
  };
  imat res = matmul(m1, m2);

  EXPECT_EQ(58, res(0, 0));
  EXPECT_EQ(64, res(0, 1));
  EXPECT_EQ(139, res(1, 0));
  EXPECT_EQ(154, res(1, 1));
}

TEST(MatrixOperationTest, FloatMatMatProd) {
  fmat m1 = {
      {1, 2, 3},
      {4, 5, 6}
  };
  fmat m2 = {
      {7, 8},
      {9, 10},
      {11, 12}
  };
  fmat res = matmul(m1, m2);

  EXPECT_EQ(58, res(0, 0));
  EXPECT_EQ(64, res(0, 1));
  EXPECT_EQ(139, res(1, 0));
  EXPECT_EQ(154, res(1, 1));
}

TEST(MatrixOperationTest, DoublePrecisionMatMatProd) {
  mat m1 = {
      {1, 2, 3},
      {4, 5, 6}
  };
  mat m2 = {
      {7, 8},
      {9, 10},
      {11, 12}
  };
  mat res = matmul(m1, m2);

  EXPECT_EQ(58, res(0, 0));
  EXPECT_EQ(64, res(0, 1));
  EXPECT_EQ(139, res(1, 0));
  EXPECT_EQ(154, res(1, 1));
}

}

#endif //MATRIX_TEST_MATRIX_OPERERATION_H
