//
// Copyright 2018-2019 The Statslabs Authors.
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

/// @file solve.h
/// @brief solve systems of linear equations.

#ifndef _SLAB_MATRIX_FNS_SOLVE_H
#define _SLAB_MATRIX_FNS_SOLVE_H

_SLAB_BEGIN_NAMESPACE

#ifndef _SLAB_USE_NO_LAPACK
template <typename T>
inline Matrix<T, 2> solve(const Matrix<T, 2> &a, const Matrix<T, 2> &b) {
  ignore(a);
  ignore(b);
  _SLAB_ERROR("solve(): unsupported element type");
}

template <>
inline Matrix<double, 2> solve(const Matrix<double, 2> &a,
                               const Matrix<double, 2> &b) {
  assert(a.n_rows() == b.n_rows());
  assert(a.n_rows() == a.n_cols());

  std::size_t n = a.n_rows();
  std::size_t nrhs = b.n_cols();
  std::size_t lda = a.n_cols();
  std::size_t ldb = b.n_cols();

  Matrix<double, 2> a_copy(a);
  Matrix<int, 1> ipiv(n);
  Matrix<double, 2> b_copy(b);

  int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, (int)n, (int)nrhs,
                           (double *)a_copy.data(), (int)lda, ipiv.data(),
                           (double *)b_copy.data(), (int)ldb);
  if (info < 0) _SLAB_ERROR("parameter i had an illegal value.");

  return b_copy;
}

template <>
inline Matrix<float, 2> solve(const Matrix<float, 2> &a,
                              const Matrix<float, 2> &b) {
  assert(a.n_rows() == b.n_rows());
  assert(a.n_rows() == a.n_cols());

  std::size_t n = a.n_rows();
  std::size_t nrhs = b.n_cols();
  std::size_t lda = a.n_cols();
  std::size_t ldb = b.n_cols();

  Matrix<float, 2> a_copy(a);
  Matrix<int, 1> ipiv(n);
  Matrix<float, 2> b_copy(b);

  int info =
      LAPACKE_sgesv(LAPACK_ROW_MAJOR, (int)n, (int)nrhs, (float *)a_copy.data(),
                    (int)lda, ipiv.data(), (float *)b_copy.data(), (int)ldb);
  if (info < 0) _SLAB_ERROR("parameter i had an illegal value.");

  return b_copy;
}

template <>
inline Matrix<std::complex<double>, 2> solve(
    const Matrix<std::complex<double>, 2> &a,
    const Matrix<std::complex<double>, 2> &b) {
  assert(a.n_rows() == b.n_rows());
  assert(a.n_rows() == a.n_cols());

  std::size_t n = a.n_rows();
  std::size_t nrhs = b.n_cols();
  std::size_t lda = a.n_cols();
  std::size_t ldb = b.n_cols();

  Matrix<std::complex<double>, 2> a_copy(a);
  Matrix<int, 1> ipiv(n);
  Matrix<std::complex<double>, 2> b_copy(b);

  int info = LAPACKE_zgesv(
      LAPACK_ROW_MAJOR, (int)n, (int)nrhs,
      reinterpret_cast<lapack_complex_double *>(a_copy.data()), (int)lda,
      ipiv.data(), reinterpret_cast<lapack_complex_double *>(b_copy.data()),
      (int)ldb);
  if (info < 0) _SLAB_ERROR("parameter i had an illegal value.");

  return b_copy;
}

template <>
inline Matrix<std::complex<float>, 2> solve(
    const Matrix<std::complex<float>, 2> &a,
    const Matrix<std::complex<float>, 2> &b) {
  assert(a.n_rows() == b.n_rows());
  assert(a.n_rows() == a.n_cols());

  std::size_t n = a.n_rows();
  std::size_t nrhs = b.n_cols();
  std::size_t lda = a.n_cols();
  std::size_t ldb = b.n_cols();

  Matrix<std::complex<float>, 2> a_copy(a);
  Matrix<int, 1> ipiv(n);
  Matrix<std::complex<float>, 2> b_copy(b);

  int info = LAPACKE_cgesv(
      LAPACK_ROW_MAJOR, (int)n, (int)nrhs,
      reinterpret_cast<lapack_complex_float *>(a_copy.data()), (int)lda,
      ipiv.data(), reinterpret_cast<lapack_complex_float *>(b_copy.data()),
      (int)ldb);
  if (info < 0) _SLAB_ERROR("parameter i had an illegal value.");

  return b_copy;
}

template <typename T>
inline Matrix<T, 2> solve(const Matrix<T, 2> &a) {
  if (a.n_rows() != a.n_cols())
    _SLAB_ERROR("solve(): matrix A should be a square matrix");

  Matrix<T, 2> b = eye<Matrix<T, 2>>(a.n_rows(), a.n_cols());
  return solve(a, b);
}
#endif

_SLAB_END_NAMESPACE

#endif  // _SLAB_MATRIX_FNS_SOLVE_H
