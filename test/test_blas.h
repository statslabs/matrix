#ifndef MATRIX_TEST_MATRIX_BLAS_H
#define MATRIX_TEST_MATRIX_BLAS_H

#include <cstring>

#include <complex>

#include "slab/matrix.h"

namespace slab {

TEST(BLASTest, LEVEL1_ASUM) {
  // value_type: double
  vec v = {6, 6, 6, 1, 2, 3};
  mat m = {{6, 1}, {6, 2}, {6, 3}};
  EXPECT_EQ(24, blas_asum(v));
  EXPECT_EQ(6, blas_asum(v.subvec(3, 5)));
  EXPECT_EQ(6, blas_asum(m.col(1)));

  // value_type: std::complex<double>
  vec cx_v_real = {6, 6, 6, 1, 2, 3};
  vec cx_v_imag = {1, 1, 1, 1, 1, 1};
  cx_vec cx_v(cx_v_real, cx_v_imag);
  EXPECT_EQ(30, blas_asum(cx_v));
}

TEST(BLASTest, LEVEL1_AXPY) {
  // value_type: double
  int a = 9;
  vec x = {6, 6, 6, 1, 2, 3};
  vec y = {6, 6, 6, 1, 2, 3};
  blas_axpy(a, x, y);

  vec y_expect = {60, 60, 60, 10, 20, 30};
  for (uword i = 0; i != y.size(); ++i) EXPECT_EQ(y_expect(i), y(i));

  vec y2 = {1, 2, 3};
  blas_axpy(a, x.subvec(3, 5), y2);

  vec y2_expect = {10, 20, 30};
  for (uword i = 0; i != y2.size(); ++i) EXPECT_EQ(y2_expect(i), y2(i));

  // value_type: std::complex<double>
  std::complex<double> cx_a(9.0, 1.0);
  vec cx_x_real = {6, 6, 6, 1, 2, 3};
  vec cx_x_imag = {1, 1, 1, 1, 1, 1};
  cx_vec cx_x(cx_x_real, cx_x_imag);
  cx_vec cx_y(cx_x_real, cx_x_imag);
  blas_axpy(cx_a, cx_x, cx_y);

  vec cx_y_real_expect = {59, 59, 59, 9, 19, 29};
  vec cx_y_imag_expect = {16, 16, 16, 11, 12, 13};
  cx_vec cx_y_expect(cx_y_real_expect, cx_y_imag_expect);
  for (uword i = 0; i != cx_y.size(); ++i) EXPECT_EQ(cx_y_expect(i), cx_y(i));
}

TEST(BLASTest, LEVEL1_COPY) {
  vec x = {6, 6, 6, 1, 2, 3}, y;

  blas_copy(x, y);

  EXPECT_EQ(6, y(0));
  EXPECT_EQ(6, y(1));
  EXPECT_EQ(6, y(2));
  EXPECT_EQ(1, y(3));
  EXPECT_EQ(2, y(4));
  EXPECT_EQ(3, y(5));

  blas_copy(x.subvec(3, 5), y);

  EXPECT_EQ(1, y(0));
  EXPECT_EQ(2, y(1));
  EXPECT_EQ(3, y(2));

  // now: y = {1, 2, 3}
  blas_copy(y, y);
  EXPECT_EQ(1, y(0));
  EXPECT_EQ(2, y(1));
  EXPECT_EQ(3, y(2));
}

TEST(BLASTest, LEVEL1_DOT) {
  vec v = {6, 6, 6, 1, 2, 3};

  EXPECT_EQ(122, blas_dot(v, v));
  EXPECT_EQ(14, blas_dot(v.subvec(3, 5), v.subvec(3, 5)));
}

TEST(BLASTest, LEVEL1_SDOT) {
  fvec v = {6, 6, 6, 1, 2, 3};

  EXPECT_EQ(123, blas_sdsdot(1.0f, v, v));
  EXPECT_EQ(122, blas_dsdot(v, v));
  EXPECT_EQ(15, blas_sdsdot(1.0f, v.subvec(3, 5), v.subvec(3, 5)));
  EXPECT_EQ(14, blas_dsdot(v.subvec(3, 5), v.subvec(3, 5)));
}

TEST(BLASTest, LEVEL1_SWAP) {
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

TEST(BLASTest, LEVEL1_SCAL) {
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

TEST(BLASTest, LEVEL1_IAMAX) {
  Matrix<double, 1> m1 = {1, 3, 2};
  Matrix<float, 1> m2 = {1, 3, 2};

  auto idx1 = blas_iamax(m1);
  auto idx2 = blas_iamax(m2);

  EXPECT_EQ(1, idx1);
  EXPECT_EQ(1, idx2);
}

TEST(BLASTest, LEVEL2_GEMV) {
  Matrix<double, 2> a1 = {{8.0, 3.0, 1.0}, {4.0, 5.0, 3.0}, {7.0, 1.0, 2.0}};

  Matrix<double, 1> x1 = {-1.0, 2.0, 1.0};
  Matrix<double, 1> y1 = {1.5, 1.0, 1.5};

  Matrix<float, 2> a2 = {{8.0, 3.0, 1.0}, {4.0, 5.0, 3.0}, {7.0, 1.0, 2.0}};

  Matrix<float, 1> x2 = {-1.0, 2.0, 1.0};
  Matrix<float, 1> y2 = {1.5, 1.0, 1.5};

  blas_gemv(CblasTrans, 1.0, a1, x1, 2.0, y1);
  blas_gemv(CblasTrans, 1.0f, a2, x2, 2.0f, y2);

  EXPECT_EQ(10, y1(0));
  EXPECT_EQ(10, y1(1));
  EXPECT_EQ(10, y1(2));
  EXPECT_EQ(10, y2(0));
  EXPECT_EQ(10, y2(1));
  EXPECT_EQ(10, y2(2));
}

TEST(BLASTest, LEVEL2_SPR) {
  Matrix<double, 1> x1 = {1.0, 2.0, 3.0};
  SymmetricMatrix<double, upper> ap1(3);

  Matrix<float, 1> x2 = {1.0, 2.0, 3.0};
  SymmetricMatrix<float, upper> ap2(3);

  for (size_t i = 0; i != 3; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      ap1(i, j) = i * 10 + j;
      ap2(i, j) = i * 10 + j;
    }
  }

  blas_spr(2.0, x1, ap1);
  blas_spr(2.0f, x2, ap2);

  EXPECT_EQ(2, ap1(0, 0));
  EXPECT_EQ(14, ap1(0, 1));
  EXPECT_EQ(26, ap1(0, 2));
  EXPECT_EQ(19, ap1(1, 1));
  EXPECT_EQ(33, ap1(1, 2));
  EXPECT_EQ(40, ap1(2, 2));
  EXPECT_EQ(2, ap2(0, 0));
  EXPECT_EQ(14, ap2(0, 1));
  EXPECT_EQ(26, ap2(0, 2));
  EXPECT_EQ(19, ap2(1, 1));
  EXPECT_EQ(33, ap2(1, 2));
  EXPECT_EQ(40, ap2(2, 2));
}

TEST(BLASTest, LEVEL2_SPR2) {
  Matrix<double, 1> x1 = {1.0, 2.0, 3.0};
  Matrix<double, 1> y1 = {1.0, 2.0, 3.0};
  SymmetricMatrix<double, upper> ap1(3);

  Matrix<float, 1> x2 = {1.0, 2.0, 3.0};
  Matrix<float, 1> y2 = {1.0, 2.0, 3.0};
  SymmetricMatrix<float, upper> ap2(3);

  for (size_t i = 0; i != 3; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      ap1(i, j) = i * 10 + j;
      ap2(i, j) = i * 10 + j;
    }
  }

  blas_spr2(1.0, x1, y1, ap1);
  blas_spr2(1.0f, x2, y2, ap2);

  EXPECT_EQ(2, ap1(0, 0));
  EXPECT_EQ(14, ap1(0, 1));
  EXPECT_EQ(26, ap1(0, 2));
  EXPECT_EQ(19, ap1(1, 1));
  EXPECT_EQ(33, ap1(1, 2));
  EXPECT_EQ(40, ap1(2, 2));
  EXPECT_EQ(2, ap2(0, 0));
  EXPECT_EQ(14, ap2(0, 1));
  EXPECT_EQ(26, ap2(0, 2));
  EXPECT_EQ(19, ap2(1, 1));
  EXPECT_EQ(33, ap2(1, 2));
  EXPECT_EQ(40, ap2(2, 2));
}

TEST(BLASTest, LEVEL3_GEMM) {
  Matrix<double, 2> a1 = {{1.0, -3.0}, {2.0, 4.0}, {1.0, -1.0}};
  Matrix<float, 2> a2 = {{1.0, -3.0}, {2.0, 4.0}, {1.0, -1.0}};

  Matrix<double, 2> b1 = a1;
  Matrix<float, 2> b2 = a2;

  Matrix<double, 2> c1 = {{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}};
  Matrix<float, 2> c2 = {{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}};

  blas_gemm(CblasNoTrans, CblasTrans, 1.0, a1, b1, 2.0, c1);
  blas_gemm(CblasNoTrans, CblasTrans, 1.0f, a2, b2, 2.0f, c2);

  EXPECT_EQ(11.0, c1(0, 0));
  EXPECT_EQ(-9.0, c1(0, 1));
  EXPECT_EQ(5.0, c1(0, 2));
  EXPECT_EQ(-9.0, c1(1, 0));
  EXPECT_EQ(21.0, c1(1, 1));
  EXPECT_EQ(-1.0, c1(1, 2));
  EXPECT_EQ(5.0, c1(2, 0));
  EXPECT_EQ(-1.0, c1(2, 1));
  EXPECT_EQ(3.0, c1(2, 2));

  EXPECT_EQ(11.0, c2(0, 0));
  EXPECT_EQ(-9.0, c2(0, 1));
  EXPECT_EQ(5.0, c2(0, 2));
  EXPECT_EQ(-9.0, c2(1, 0));
  EXPECT_EQ(21.0, c2(1, 1));
  EXPECT_EQ(-1.0, c2(1, 2));
  EXPECT_EQ(5.0, c2(2, 0));
  EXPECT_EQ(-1.0, c2(2, 1));
  EXPECT_EQ(3.0, c2(2, 2));
}

}  // namespace slab

#endif  // MATRIX_TEST_MATRIX_BLAS_H
