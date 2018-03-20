//
// Created by Yi Pan (Institute of Cancer and Genomic Sciences) on 19/03/2018.
//

#ifndef MATRIX_TEST_MATRIX_BLAS_H
#define MATRIX_TEST_MATRIX_BLAS_H

#include "slab/matrix.h"

namespace slab {

TEST(BLASlevel1Test, SWAP) {
  Matrix<double, 1> m1 = {4, 5, 6};
  Matrix<double, 1> m1s = {1, 2, 3};

  Matrix<float, 1> m2 = {4, 5, 6};
  Matrix<float, 1> m2s = {1, 2, 3};

  blas_swap(m1, m1s);
  blas_swap(m2, m2s);

  EXPECT_EQ(1, m1(0));
  EXPECT_EQ(2, m1(1));
  EXPECT_EQ(3, m1(2));
  EXPECT_EQ(4, m1s(0));
  EXPECT_EQ(5, m1s(1));
  EXPECT_EQ(6, m1s(2));

  EXPECT_EQ(1, m2(0));
  EXPECT_EQ(2, m2(1));
  EXPECT_EQ(3, m2(2));
  EXPECT_EQ(4, m2s(0));
  EXPECT_EQ(5, m2s(1));
  EXPECT_EQ(6, m2s(2));
}

TEST(BLASlevel1Test, SCAL) {
  double a1 = 0.1;
  Matrix<double, 1> m1 = {10, 20, 30};

  float a2 = 0.1;
  Matrix<float, 1> m2 = {10, 20, 30};

  blas_scal(a1, m1);
  blas_scal(a2, m2);

  EXPECT_EQ(1, m1(0));
  EXPECT_EQ(2, m1(1));
  EXPECT_EQ(3, m1(2));

  EXPECT_EQ(1, m2(0));
  EXPECT_EQ(2, m2(1));
  EXPECT_EQ(3, m2(2));
}

TEST(BLASlevel1Test, COPY) {
  Matrix<double, 1> m1 = {1, 2, 3}, m1c;
  Matrix<float, 1> m2 = {1, 2, 3}, m2c;

  blas_copy(m1, m1c);
  blas_copy(m2, m2c);

  EXPECT_EQ(1, m1c(0));
  EXPECT_EQ(2, m1c(1));
  EXPECT_EQ(3, m1c(2));

  EXPECT_EQ(1, m2c(0));
  EXPECT_EQ(2, m2c(1));
  EXPECT_EQ(3, m2c(2));
}

TEST(BLASlevel1Test, AXPY) {
  double a1 = 9.0;
  Matrix<double, 1> x1 = {1, 2, 3};
  Matrix<double, 1> y1 = {1, 2, 3};

  float a2 = 9.0;
  Matrix<float, 1> x2 = {1, 2, 3};
  Matrix<float, 1> y2 = {1, 2, 3};

  blas_axpy(a1, x1, y1);
  blas_axpy(a2, x2, y2);

  EXPECT_EQ(10, y1(0));
  EXPECT_EQ(20, y1(1));
  EXPECT_EQ(30, y1(2));
  EXPECT_EQ(10, y2(0));
  EXPECT_EQ(20, y2(1));
  EXPECT_EQ(30, y2(2));
}

TEST(BLASlevel1Test, DOT) {
  Matrix<double, 1> m1 = {1, 2, 3};
  Matrix<float, 1> m2 = {1, 2, 3};

  EXPECT_EQ(14, blas_dot(m1, m1));
  EXPECT_EQ(14, blas_dot(m2, m2));
}

TEST(BLASlevel1Test, IAMAX) {
  Matrix<double, 1> m1 = {1, 3, 2};
  Matrix<float, 1> m2 = {1, 3, 2};

  auto idx1 = blas_iamax(m1);
  auto idx2 = blas_iamax(m2);

  EXPECT_EQ(1, idx1);
  EXPECT_EQ(1, idx2);
}

}

#endif //MATRIX_TEST_MATRIX_BLAS_H
