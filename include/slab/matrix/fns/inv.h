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

/// @file inv.h
/// @brief inverse of general square matrix.

#ifndef SLAB_MATRIX_FNS_INV_H_
#define SLAB_MATRIX_FNS_INV_H_

namespace slab {

template <typename T>
inline Matrix<T, 2> inverse(const Matrix<T, 2> &a) {
  assert(a.n_rows() == a.n_cols());

  int info;
  const int m = a.n_rows();
  const int n = a.n_cols();
  const int lda = a.n_cols();
  Matrix<int, 1> ipiv(n);

  Matrix<T, 2> res = a;
  if (is_double<T>::value) {
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, (double *)res.data(), lda,
                          ipiv.data());
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, (double *)res.data(), lda, ipiv.data());
  } else if (is_float<T>::value) {
    info = LAPACKE_sgetrf(LAPACK_ROW_MAJOR, m, n, (float *)res.data(), lda,
                          ipiv.data());
    LAPACKE_sgetri(LAPACK_ROW_MAJOR, n, (float *)res.data(), lda, ipiv.data());
  } else if (is_complex_double<T>::value) {
    info =
        LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m, n,
                       (lapack_complex_double *)res.data(), lda, ipiv.data());
    LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, (lapack_complex_double *)res.data(),
                   lda, ipiv.data());
  } else if (is_complex_float<T>::value) {
    info = LAPACKE_cgetrf(LAPACK_ROW_MAJOR, m, n,
                          (lapack_complex_float *)res.data(), lda, ipiv.data());
    LAPACKE_cgetri(LAPACK_ROW_MAJOR, n, (lapack_complex_float *)res.data(), lda,
                   ipiv.data());
  } else {
    err_msg("inverse(): unspported element type.");
  }

  return res;
}

template <typename T>
inline bool inv(Matrix<T, 2> &b, const Matrix<T, 2> &a) {
  return true;
}

template <>
inline bool inv(Matrix<double, 2> &b, const Matrix<double, 2> &a) {
  b = eye<Matrix<double, 2>>(a.n_rows(), a.n_cols());

  std::size_t n = a.n_rows();
  std::size_t nrhs = b.n_cols();
  std::size_t lda = a.n_cols();
  std::size_t ldb = b.n_cols();

  Matrix<double, 2> a_copy = a;
  Matrix<int, 1> ipiv(n);

  int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, (int)n, (int)nrhs,
                           (double *)a_copy.data(), (int)lda, ipiv.data(),
                           (double *)b.data(), (int)ldb);
  if (info) return false;

  return true;
}

}  // namespace slab

#endif  // SLAB_MATRIX_FNS_INV_H_
