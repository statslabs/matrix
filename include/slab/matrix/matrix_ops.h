//
// Copyright 2018 The Statslabs Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// -----------------------------------------------------------------------------
// operations.h
// -----------------------------------------------------------------------------
//
#ifndef SLAB_MATRIX_OPERATIONS_H_
#define SLAB_MATRIX_OPERATIONS_H_

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <complex>

#ifdef USE_MKL
#include "mkl.h"
#else
extern "C" {
#include "cblas.h"
}
#endif

#include "slab/matrix/error.h"
#include "slab/matrix/matrix.h"
#include "slab/matrix/matrix_base.h"
#include "slab/matrix/matrix_ref.h"
#include "slab/matrix/support.h"

namespace slab {

template <typename M1, typename M2>
inline Enable_if<Matrix_type<M1>() && Matrix_type<M2>(), bool> operator==(
    const M1 &a, const M2 &b) {
  assert(same_extents(a.descriptor(), b.descriptor()));
  return std::equal(a.begin(), a.end(), b.begin());
}

template <typename M1, typename M2>
inline Enable_if<Matrix_type<M1>() && Matrix_type<M2>(), bool> operator!=(
    const M1 &a, const M2 &b) {
  return !(a == b);
}

// Scalar Addtion
//
// res = X + val or res = val + X

template <typename T, std::size_t N>
Matrix<T, N> operator+(const Matrix<T, N> &x, const T &val) {
  Matrix<T, N> res = x;
  res += val;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator+(const MatrixRef<T, N> &x, const T &val) {
  Matrix<T, N> res = x;
  res += val;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator+(const T &val, const Matrix<T, N> &x) {
  Matrix<T, N> res = x;
  res += val;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator+(const T &val, const MatrixRef<T, N> &x) {
  Matrix<T, N> res = x;
  res += val;
  return res;
}

// Scalar Subtraction
//
// res = X - val

template <typename T, std::size_t N>
Matrix<T, N> operator-(const Matrix<T, N> &x, const T &val) {
  Matrix<T, N> res = x;
  res -= val;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator-(const MatrixRef<T, N> &x, const T &val) {
  Matrix<T, N> res = x;
  res -= val;
  return res;
}

// Scalar Multiplication
//
// res = X * val or res = val * X

template <typename T, std::size_t N>
Matrix<T, N> operator*(const Matrix<T, N> &x, const T &val) {
  Matrix<T, N> res = x;
  res *= val;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator*(const MatrixRef<T, N> &x, const T &val) {
  Matrix<T, N> res = x;
  res *= val;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator*(const T &val, const Matrix<T, N> &x) {
  Matrix<T, N> res = x;
  res *= val;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator*(const T &val, const MatrixRef<T, N> &x) {
  Matrix<T, N> res = x;
  res *= val;
  return res;
}

// Scalar Division
//
// res = X / val

template <typename T, std::size_t N>
Matrix<T, N> operator/(const Matrix<T, N> &x, const T &val) {
  Matrix<T, N> res = x;
  res /= val;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator/(const MatrixRef<T, N> &x, const T &val) {
  Matrix<T, N> res = x;
  res /= val;
  return res;
}

// Scalar Modulus
//
// res = X % val

template <typename T, std::size_t N>
Matrix<T, N> operator%(const Matrix<T, N> &x, const T &val) {
  Matrix<T, N> res = x;
  res %= val;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator%(const MatrixRef<T, N> &x, const T &val) {
  Matrix<T, N> res = x;
  res %= val;
  return res;
}

// Matrix Addtion
//
// res = A + B

template <typename T, std::size_t N>
Matrix<T, N> operator+(const Matrix<T, N> &a, const Matrix<T, N> &b) {
  Matrix<T, N> res = a;
  res += b;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator+(const MatrixRef<T, N> &a, const MatrixRef<T, N> &b) {
  Matrix<T, N> res = a;
  res += b;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator+(const Matrix<T, N> &a, const MatrixRef<T, N> &b) {
  Matrix<T, N> res = a;
  res += b;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator+(const MatrixRef<T, N> &a, const Matrix<T, N> &b) {
  Matrix<T, N> res = a;
  res += b;
  return res;
}

// Matrix Subtraction
//
// res = A - B

template <typename T, std::size_t N>
Matrix<T, N> operator-(const Matrix<T, N> &a, const Matrix<T, N> &b) {
  Matrix<T, N> res = a;
  res -= b;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator-(const MatrixRef<T, N> &a, const MatrixRef<T, N> &b) {
  Matrix<T, N> res = a;
  res -= b;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator-(const Matrix<T, N> &a, const MatrixRef<T, N> &b) {
  Matrix<T, N> res = a;
  res -= b;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator-(const MatrixRef<T, N> &a, const Matrix<T, N> &b) {
  Matrix<T, N> res = a;
  res -= b;
  return res;
}

// Element-wise Multiplication
//
// res = A * B

template <typename T, std::size_t N>
Matrix<T, N> operator*(const Matrix<T, N> &a, const Matrix<T, N> &b) {
  Matrix<T, N> res = a;
  res *= b;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator*(const MatrixRef<T, N> &a, const MatrixRef<T, N> &b) {
  Matrix<T, N> res = a;
  res *= b;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator*(const Matrix<T, N> &a, const MatrixRef<T, N> &b) {
  Matrix<T, N> res = a;
  res *= b;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator*(const MatrixRef<T, N> &a, const Matrix<T, N> &b) {
  Matrix<T, N> res = a;
  res *= b;
  return res;
}

// Element-wise Division
//
// res = A / B

template <typename T, std::size_t N>
Matrix<T, N> operator/(const Matrix<T, N> &a, const Matrix<T, N> &b) {
  Matrix<T, N> res = a;
  res /= b;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator/(const MatrixRef<T, N> &a, const MatrixRef<T, N> &b) {
  Matrix<T, N> res = a;
  res /= b;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator/(const Matrix<T, N> &a, const MatrixRef<T, N> &b) {
  Matrix<T, N> res = a;
  res /= b;
  return res;
}

template <typename T, std::size_t N>
Matrix<T, N> operator/(const MatrixRef<T, N> &a, const Matrix<T, N> &b) {
  Matrix<T, N> res = a;
  res /= b;
  return res;
}

//
// Matrix Multiplication
//

template <typename T>
inline Matrix<T, 2> matmul(const MatrixBase<T, 1> &a,
                           const MatrixBase<T, 2> &b) {
  assert(b.n_rows() == 1);

  Matrix<T, 2> mat_a(a.n_rows(), 1);
  for (std::size_t i = 0; i != a.n_rows(); ++i) mat_a(i, 0) = a(i);

  return matmul(mat_a, b);
}

template <typename T>
inline Matrix<T, 1> matmul(const MatrixBase<T, 2> &a,
                           const MatrixBase<T, 1> &x) {
  assert(a.extent(1) == x.extent(0));

  const std::size_t m = a.n_rows();
  const std::size_t n = a.n_cols();
  Matrix<T, 1> y(m);

  for (std::size_t i = 0; i != m; ++i)
    for (std::size_t j = 0; j != n; ++j) y(i) += a(i, j) * x(j);

  return y;
}

template <>
inline Matrix<double, 1> matmul(const MatrixBase<double, 2> &a,
                                const MatrixBase<double, 1> &x) {
  assert(a.extent(1) == x.extent(0));
  const std::size_t m = a.n_rows();
  const std::size_t n = a.n_cols();
  const std::size_t lda = n;
  const std::size_t incx = x.descriptor().strides[0];
  const std::size_t incy = 1;

  Matrix<double, 1> y(m);
  cblas_dgemv(CblasRowMajor, CblasNoTrans, (const int) m, (const int) n,
	      (const double)1.0,
              (const double *)(a.data() + a.descriptor().start),
	      (const int) lda,
              (const double *)(x.data() + x.descriptor().start),
	      (const int) incx,
              (const double)0.0, (double *)y.data(),
	      (const int) incy);

  return y;
}

template <>
inline Matrix<float, 1> matmul(const MatrixBase<float, 2> &a,
                               const MatrixBase<float, 1> &x) {
  assert(a.extent(1) == x.extent(0));
  const std::size_t m = a.n_rows();
  const std::size_t n = a.n_cols();
  const std::size_t lda = n;
  const std::size_t incx = x.descriptor().strides[0];
  const std::size_t incy = 1;

  Matrix<float, 1> y(m);
  cblas_sgemv(CblasRowMajor, CblasNoTrans, (const int) m, (const int) n,
	      (const float)1.0,
              (const float *)(a.data() + a.descriptor().start),
	      (const int) lda,
              (const float *)(x.data() + x.descriptor().start),
	      (const int) incx,
              (const float)0.0, (float *)y.data(),
	      (const int) incy);

  return y;
}

template <typename T>
inline Matrix<T, 2> matmul(const MatrixBase<T, 2> &a,
                           const MatrixBase<T, 2> &b) {
  assert(a.extent(1) == b.extent(0));

  const std::size_t m = a.n_rows();
  const std::size_t n = b.n_cols();
  const std::size_t k = a.n_cols();

  Matrix<T, 2> c(m, n);

  for (std::size_t i = 0; i != m; ++i) {
    for (std::size_t j = 0; j != n; ++j) {
      for (std::size_t idx = 0; idx != k; ++idx) {
        c(i, j) += a(i, idx) * b(idx, j);
      }
    }
  }

  return c;
}

template <>
inline Matrix<double, 2> matmul(const MatrixBase<double, 2> &a,
                                const MatrixBase<double, 2> &b) {
  assert(a.extent(1) == b.extent(0));

  const std::size_t m = a.n_rows();
  const std::size_t n = b.n_cols();
  const std::size_t k = a.n_cols();

  const std::size_t lda = a.n_cols();
  const std::size_t ldb = b.n_cols();
  const std::size_t ldc = b.n_cols();

  Matrix<double, 2> c(m, n);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	      (const int) m, (const int) n, (const int) k,
              (const double)1.0,
              (const double *)(a.data() + a.descriptor().start),
	      (const int) lda,
              (const double *)(b.data() + b.descriptor().start),
	      (const int) ldb,
              (const double)0.0, (double *)c.data(),
	      (const int) ldc);

  return c;
}

template <>
inline Matrix<float, 2> matmul(const MatrixBase<float, 2> &a,
                               const MatrixBase<float, 2> &b) {
  assert(a.extent(1) == b.extent(0));

  const std::size_t m = a.n_rows();
  const std::size_t n = b.n_cols();
  const std::size_t k = a.n_cols();

  const std::size_t lda = a.n_cols();
  const std::size_t ldb = b.n_cols();
  const std::size_t ldc = b.n_cols();

  Matrix<float, 2> c(m, n);
  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	      (const int) m, (const int) n, (const int) k,
              (const float)1.0,
              (const float *)(a.data() + a.descriptor().start),
	      (const int) lda,
              (const float *)(b.data() + b.descriptor().start),
	      (const int) ldb,
              (const float)0.0, (float *)c.data(),
	      (const int) ldc);

  return c;
}

template <typename T>
inline const Matrix<T, 1> &matmul_n(const Matrix<T, 1> &x) {
  return x;
}

template <typename T>
inline const Matrix<T, 2> &matmul_n(const Matrix<T, 2> &x) {
  return x;
}

template <typename T, typename... Args>
inline auto matmul_n(const Matrix<T, 2> &x, Args... args)
    -> decltype(matmul(x, matmul_n(args...))) {
  return matmul(x, matmul_n(args...));
}

template <typename T>
inline Matrix<T, 2> diag(const Matrix<T, 1> &x) {
  Matrix<T, 2> res(x.size(), x.size());
  res.diag() = x;

  return res;
}

template <typename T>
inline T dot(const MatrixBase<T, 1> &a, const MatrixBase<T, 1> &b) {
  assert(a.size() == b.size());

  T res = T{0};
  for (std::size_t idx = 0; idx != a.size(); ++idx) {
    res += a(idx) * b(idx);
  }

  return res;
}

}  // namespace slab

#endif  // SLAB_MATRIX_OPERATIONS_H_
