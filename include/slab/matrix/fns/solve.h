//
// Copyright 2019 The Statslabs Authors.
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

#ifndef SLAB_MATRIX_FNS_SOLVE_H_
#define SLAB_MATRIX_FNS_SOLVE_H_

namespace slab {

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

}  // namespace slab

#endif  // SLAB_MATRIX_FNS_SOLVE_H_
