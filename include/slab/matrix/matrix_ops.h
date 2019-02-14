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

#include "slab/matrix/blas_interface.h"
#include "slab/matrix/lapack_interface.h"

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
  const int m = a.n_rows();
  const int n = a.n_cols();
  const int lda = n;
  const int incx = x.descriptor().strides[0];
  const int incy = 1;

  Matrix<double, 1> y(m);
  cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, (const double)1.0,
              (const double *)(a.data() + a.descriptor().start), lda,
              (const double *)(x.data() + x.descriptor().start), incx,
              (const double)0.0, (double *)y.data(), incy);

  return y;
}

template <>
inline Matrix<float, 1> matmul(const MatrixBase<float, 2> &a,
                               const MatrixBase<float, 1> &x) {
  assert(a.extent(1) == x.extent(0));
  const int m = a.n_rows();
  const int n = a.n_cols();
  const int lda = n;
  const int incx = x.descriptor().strides[0];
  const int incy = 1;

  Matrix<float, 1> y(m);
  cblas_sgemv(CblasRowMajor, CblasNoTrans, m, n, (const float)1.0,
              (const float *)(a.data() + a.descriptor().start), lda,
              (const float *)(x.data() + x.descriptor().start), incx,
              (const float)0.0, (float *)y.data(), incy);

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

  const int m = a.n_rows();
  const int n = b.n_cols();
  const int k = a.n_cols();

  const int lda = a.n_cols();
  const int ldb = b.n_cols();
  const int ldc = b.n_cols();

  Matrix<double, 2> c(m, n);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              (const double)1.0,
              (const double *)(a.data() + a.descriptor().start), lda,
              (const double *)(b.data() + b.descriptor().start), ldb,
              (const double)0.0, (double *)c.data(), ldc);

  return c;
}

template <>
inline Matrix<float, 2> matmul(const MatrixBase<float, 2> &a,
                               const MatrixBase<float, 2> &b) {
  assert(a.extent(1) == b.extent(0));

  const int m = a.n_rows();
  const int n = b.n_cols();
  const int k = a.n_cols();

  const int lda = a.n_cols();
  const int ldb = b.n_cols();
  const int ldc = b.n_cols();

  Matrix<float, 2> c(m, n);
  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              (const float)1.0,
              (const float *)(a.data() + a.descriptor().start), lda,
              (const float *)(b.data() + b.descriptor().start), ldb,
              (const float)0.0, (float *)c.data(), ldc);

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

template <typename T, std::size_t N, typename... Args>
inline auto reshape(const Matrix<T, N> &x, Args... args)
    -> decltype(Matrix<T, sizeof...(args)>()) {
  Matrix<T, sizeof...(args)> res(args...);
  std::copy(x.begin(), x.end(), res.begin());

  return res;
}

template <typename T, std::size_t N>
inline Matrix<T, 1> vectorise(const Matrix<T, N> &x) {
  return reshape(x, x.size());
}

// join_vecs()
template <typename T>
inline Matrix<T, 1> join_vecs(const Matrix<T, 1> &a, const Matrix<T, 1> &b) {
  Matrix<T, 1> res(a.n_rows() + b.n_rows());
  res(slice{0, a.n_rows()}) = a;
  res(slice{a.n_rows(), b.n_rows()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 1> join_vecs(const MatrixRef<T, 1> &a,
                              const MatrixRef<T, 1> &b) {
  Matrix<T, 1> res(a.n_rows() + b.n_rows());
  res(slice{0, a.n_rows()}) = a;
  res(slice{a.n_rows(), b.n_rows()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 1> join_vecs(const Matrix<T, 1> &a, const MatrixRef<T, 1> &b) {
  Matrix<T, 1> res(a.n_rows() + b.n_rows());
  res(slice{0, a.n_rows()}) = a;
  res(slice{a.n_rows(), b.n_rows()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 1> join_vecs(const MatrixRef<T, 1> &a, const Matrix<T, 1> &b) {
  Matrix<T, 1> res(a.n_rows() + b.n_rows());
  res(slice{0, a.n_rows()}) = a;
  res(slice{a.n_rows(), b.n_rows()}) = b;

  return res;
}

template <typename T, typename... Args>
inline Matrix<T, 1> join_vecs(const Matrix<T, 1> &x, Args... args) {
  return join_vecs(x, join_vecs(args...));
}

template <typename T, typename... Args>
inline Matrix<T, 1> join_vecs(const MatrixRef<T, 1> &x, Args... args) {
  return join_vecs(x, join_vecs(args...));
}

// join_rows()

template <typename T>
inline Matrix<T, 2> join_rows(const Matrix<T, 2> &a, const Matrix<T, 2> &b) {
  Matrix<T, 2> res;
  if (a.empty() && b.empty())
    return res;
  else if (a.empty())
    return b;
  else if (b.empty())
    return a;
  else if (a.n_rows() != b.n_rows())
    err_quit("joint_rows(): inconsistent number of rows");

  res = slab::zeros<Matrix<T, 2>>(a.n_rows(), a.n_cols() + b.n_cols());
  res(slice{0, a.n_rows()}, slice{0, a.n_cols()}) = a;
  res(slice{0, a.n_rows()}, slice{a.n_cols(), b.n_cols()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 2> join_rows(const MatrixRef<T, 2> &a,
                              const MatrixRef<T, 2> &b) {
  assert(a.n_rows() == b.n_rows());

  Matrix<T, 2> res(a.n_rows(), a.n_cols() + b.n_cols());
  res(slice{0, a.n_rows()}, slice{0, a.n_cols()}) = a;
  res(slice{0, a.n_rows()}, slice{a.n_cols(), b.n_cols()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 2> join_rows(const Matrix<T, 2> &a, const MatrixRef<T, 2> &b) {
  assert(a.n_rows() == b.n_rows());

  Matrix<T, 2> res(a.n_rows(), a.n_cols() + b.n_cols());
  res(slice{0, a.n_rows()}, slice{0, a.n_cols()}) = a;
  res(slice{0, a.n_rows()}, slice{a.n_cols(), b.n_cols()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 2> join_rows(const MatrixRef<T, 2> &a, const Matrix<T, 2> &b) {
  assert(a.n_rows() == b.n_rows());

  Matrix<T, 2> res(a.n_rows(), a.n_cols() + b.n_cols());
  res(slice{0, a.n_rows()}, slice{0, a.n_cols()}) = a;
  res(slice{0, a.n_rows()}, slice{a.n_cols(), b.n_cols()}) = b;

  return res;
}

// join_cols()

template <typename T>
inline Matrix<T, 2> join_cols(const Matrix<T, 2> &a, const Matrix<T, 2> &b) {
  Matrix<T, 2> res;
  if (a.empty() && b.empty())
    return res;
  else if (a.empty())
    return b;
  else if (b.empty())
    return a;
  else if (a.n_cols() != b.n_cols())
    err_quit("joint_rows(): inconsistent number of columns");

  res = slab::zeros<Matrix<T, 2>>(a.n_rows() + b.n_rows(), a.n_cols());
  res(slice{0, a.n_rows()}, slice{0, a.n_cols()}) = a;
  res(slice{a.n_rows(), b.n_rows()}, slice{0, a.n_cols()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 2> join_cols(const MatrixRef<T, 2> &a,
                              const MatrixRef<T, 2> &b) {
  assert(a.n_cols() == b.n_cols());

  Matrix<T, 2> res(a.n_rows() + b.n_rows(), a.n_cols());
  res(slice{0, a.n_rows()}, slice{0, a.n_cols()}) = a;
  res(slice{a.n_rows(), b.n_rows()}, slice{0, a.n_cols()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 2> join_cols(const Matrix<T, 2> &a, const MatrixRef<T, 2> &b) {
  assert(a.n_cols() == b.n_cols());

  Matrix<T, 2> res(a.n_rows() + b.n_rows(), a.n_cols());
  res(slice{0, a.n_rows()}, slice{0, a.n_cols()}) = a;
  res(slice{a.n_rows(), b.n_rows()}, slice{0, a.n_cols()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 2> join_cols(const MatrixRef<T, 2> &a, const Matrix<T, 2> &b) {
  assert(a.n_cols() == b.n_cols());

  Matrix<T, 2> res(a.n_rows() + b.n_rows(), a.n_cols());
  res(slice{0, a.n_rows()}, slice{0, a.n_cols()}) = a;
  res(slice{a.n_rows(), b.n_rows()}, slice{0, a.n_cols()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 2> kron(const Matrix<T, 2> &a, const Matrix<T, 2> &b) {
  const std::size_t a_rows = a.n_rows();
  const std::size_t a_cols = a.n_cols();
  const std::size_t b_rows = b.n_rows();
  const std::size_t b_cols = b.n_cols();

  Matrix<T, 2> res(a_rows * b_rows, a_cols * b_cols);
  for (std::size_t j = 0; j != a_cols; ++j) {
    for (std::size_t i = 0; i != a_rows; ++i) {
      res(slice{i * b_rows, b_rows}, slice{j * b_cols, b_cols}) = a(i, j) * b;
    }
  }

  return res;
}

template <typename T>
inline Matrix<T, 2> solve(const Matrix<T, 2> &a, const Matrix<T, 2> &b) {
  err_quit("solve(): unsupported element type");
}

template <>
inline Matrix<double, 2> solve(const Matrix<double, 2> &a,
                               const Matrix<double, 2> &b) {
  assert(a.n_rows() == b.n_rows());
  assert(a.n_rows() == a.n_cols());

  int n = a.n_rows();
  int nrhs = b.n_cols();
  int lda = a.n_cols();
  int ldb = b.n_cols();

  Matrix<double, 2> a_copy(a);
  Matrix<int, 1> ipiv(n);
  Matrix<double, 2> b_copy(b);

  int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, (double *)a_copy.data(),
                           lda, ipiv.data(), (double *)b_copy.data(), ldb);

  return b_copy;
}

template <>
inline Matrix<float, 2> solve(const Matrix<float, 2> &a,
                              const Matrix<float, 2> &b) {
  assert(a.n_rows() == b.n_rows());
  assert(a.n_rows() == a.n_cols());

  int n = a.n_rows();
  int nrhs = b.n_cols();
  int lda = a.n_cols();
  int ldb = b.n_cols();

  Matrix<float, 2> a_copy(a);
  Matrix<int, 1> ipiv(n);
  Matrix<float, 2> b_copy(b);

  int info = LAPACKE_sgesv(LAPACK_ROW_MAJOR, n, nrhs, (float *)a_copy.data(),
                           lda, ipiv.data(), (float *)b_copy.data(), ldb);

  return b_copy;
}

template <>
inline Matrix<std::complex<double>, 2> solve(
    const Matrix<std::complex<double>, 2> &a,
    const Matrix<std::complex<double>, 2> &b) {
  assert(a.n_rows() == b.n_rows());
  assert(a.n_rows() == a.n_cols());

  int n = a.n_rows();
  int nrhs = b.n_cols();
  int lda = a.n_cols();
  int ldb = b.n_cols();

  Matrix<std::complex<double>, 2> a_copy(a);
  Matrix<int, 1> ipiv(n);
  Matrix<std::complex<double>, 2> b_copy(b);

  int info = LAPACKE_zgesv(
      LAPACK_ROW_MAJOR, n, nrhs,
      reinterpret_cast<lapack_complex_double *>(a_copy.data()), lda,
      ipiv.data(), reinterpret_cast<lapack_complex_double *>(b_copy.data()),
      ldb);

  return b_copy;
}

template <>
inline Matrix<std::complex<float>, 2> solve(
    const Matrix<std::complex<float>, 2> &a,
    const Matrix<std::complex<float>, 2> &b) {
  assert(a.n_rows() == b.n_rows());
  assert(a.n_rows() == a.n_cols());

  int n = a.n_rows();
  int nrhs = b.n_cols();
  int lda = a.n_cols();
  int ldb = b.n_cols();

  Matrix<std::complex<float>, 2> a_copy(a);
  Matrix<int, 1> ipiv(n);
  Matrix<std::complex<float>, 2> b_copy(b);

  int info = LAPACKE_cgesv(
      LAPACK_ROW_MAJOR, n, nrhs,
      reinterpret_cast<lapack_complex_float *>(a_copy.data()), lda, ipiv.data(),
      reinterpret_cast<lapack_complex_float *>(b_copy.data()), ldb);

  return b_copy;
}

template <typename T>
inline Matrix<T, 2> solve(const Matrix<T, 2> &a) {
  if (a.n_rows() != a.n_cols())
    err_quit("solve(): matrix A should be a square matrix");

  Matrix<T, 2> b = eye<Matrix<T, 2>>(a.n_rows(), a.n_cols());
  return solve(a, b);
}

template <typename T>
inline bool inv(Matrix<T, 2> &b, const Matrix<T, 2> &a) {
  return true;
}

template <>
inline bool inv(Matrix<double, 2> &b, const Matrix<double, 2> &a) {
  b = eye<Matrix<double, 2>>(a.n_rows(), a.n_cols());

  int n = a.n_rows();
  int nrhs = b.n_cols();
  int lda = a.n_cols();
  int ldb = b.n_cols();

  Matrix<double, 2> a_copy = a;
  Matrix<int, 1> ipiv(n);

  int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, (double *)a_copy.data(),
                           lda, ipiv.data(), (double *)b.data(), ldb);
  if (info) return false;

  return true;
}

template <typename T>
inline bool pinv(Matrix<T, 2> &a_inv, const Matrix<T, 2> &a) {
  int m = a.n_rows();
  int n = a.n_cols();
  int k = std::min(m, n);

  Matrix<T, 2> a_copy = a, u, vt;
  Matrix<T, 1> s, superb;

  int info = lapack_gesvd('S', 'S', a_copy, s, u, vt, superb);
  if (info) return false;

  for (int i = 0; i < k; i++) {
    double ss;
    if (s[i] > 1.0e-9)
      ss = 1.0 / s[i];
    else
      ss = s[i];

    Matrix<T, 1> ui = u.col(i);
    blas_scal(ss, ui);
    u.col(i) = ui;
  }

  a_inv = zeros<Matrix<T, 2>>(n, m);
  T alpha = 1.0;
  T beta = 0.0;
  blas_gemm(CblasTrans, CblasTrans, alpha, vt, u, beta, a_inv);

  return true;
}

template <typename T, std::size_t N>
inline Matrix<T, N> exp(const Matrix<T, N> &x) {
  Matrix<T, N> res = x;
  res.apply([](T &a) { a = std::exp(a); });

  return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> exp(const MatrixRef<T, N> &x) {
  Matrix<T, N> res = x;
  res.apply([](T &a) { a = std::exp(a); });

  return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> log(const Matrix<T, N> &x) {
  Matrix<T, N> res = x;
  res.apply([](T &a) { a = std::log(a); });

  return res;
}

template <typename T, std::size_t N>
inline Matrix<T, N> log(const MatrixRef<T, N> &x) {
  Matrix<T, N> res = x;
  res.apply([](T &a) { a = std::log(a); });

  return res;
}

template <typename T, typename T1, std::size_t N>
inline Matrix<T, N> pow(const Matrix<T, N> &x, const T1 &val) {
  static_assert(Convertible<T1, T>(), "pow(): incompatible element types");

  Matrix<T, N> res = x;
  res.apply([&](T &a) { a = std::pow(a, static_cast<T>(val)); });

  return res;
}

template <typename T, typename T1, std::size_t N>
inline Matrix<T, N> pow(const MatrixRef<T, N> &x, const T1 &val) {
  static_assert(Convertible<T1, T>(), "pow(): incompatible element types");

  Matrix<T, N> res = x;
  res.apply([&](T &a) { a = std::pow(a, static_cast<T>(val)); });

  return res;
}

template <typename T>
inline Matrix<T, 2> chol(const Matrix<T, 2> &x) {
  assert(x.n_rows() == x.n_cols());

  Matrix<T, 2> x_copy(x);
  int info = lapack_potrf(x_copy);
  if (info) err_quit("chol(): unsuccessful");

  Matrix<T, 2> res = zeros<Matrix<T, 2>>(x.n_rows(), x.n_cols());
  for (std::size_t i = 0; i != x.n_rows(); ++i) {
    for (std::size_t j = i; j != x.n_cols(); ++j) {
      res(i, j) = x_copy(i, j);
    }
  }

  return res;
}

}  // namespace slab

#endif  // SLAB_MATRIX_OPERATIONS_H_
