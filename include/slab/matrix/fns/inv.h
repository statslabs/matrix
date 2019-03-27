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
inline bool inverse(Matrix<T, 2> &b, const Matrix<T, 2> &a) {
  assert(a.n_rows() == a.n_cols());

  int info;
  const int m = a.n_rows();
  const int n = a.n_cols();
  const int lda = a.n_cols();
  Matrix<int, 1> ipiv(n);

  b = a;
  if (is_double<T>::value) {
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, (double *)b.data(), lda,
                          ipiv.data());
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, (double *)b.data(), lda, ipiv.data());
  } else if (is_float<T>::value) {
    info = LAPACKE_sgetrf(LAPACK_ROW_MAJOR, m, n, (float *)b.data(), lda,
                          ipiv.data());
    LAPACKE_sgetri(LAPACK_ROW_MAJOR, n, (float *)b.data(), lda, ipiv.data());
  } else if (is_complex_double<T>::value) {
    info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m, n,
                          (lapack_complex_double *)b.data(), lda, ipiv.data());
    LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, (lapack_complex_double *)b.data(), lda,
                   ipiv.data());
  } else if (is_complex_float<T>::value) {
    info = LAPACKE_cgetrf(LAPACK_ROW_MAJOR, m, n,
                          (lapack_complex_float *)b.data(), lda, ipiv.data());
    LAPACKE_cgetri(LAPACK_ROW_MAJOR, n, (lapack_complex_float *)b.data(), lda,
                   ipiv.data());
  } else {
    err_msg("inverse(): unspported element type.");
  }

  if (info) return false;

  return true;
}

template <typename T>
inline bool inverse(Matrix<T, 2> &b, const MatrixRef<T, 2> &a) {
  assert(a.n_rows() == a.n_cols());

  int info;
  const int m = a.n_rows();
  const int n = a.n_cols();
  const int lda = a.n_cols();
  Matrix<int, 1> ipiv(n);

  b = a;
  if (is_double<T>::value) {
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, (double *)b.data(), lda,
                          ipiv.data());
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, (double *)b.data(), lda, ipiv.data());
  } else if (is_float<T>::value) {
    info = LAPACKE_sgetrf(LAPACK_ROW_MAJOR, m, n, (float *)b.data(), lda,
                          ipiv.data());
    LAPACKE_sgetri(LAPACK_ROW_MAJOR, n, (float *)b.data(), lda, ipiv.data());
  } else if (is_complex_double<T>::value) {
    info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m, n,
                          (lapack_complex_double *)b.data(), lda, ipiv.data());
    LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, (lapack_complex_double *)b.data(), lda,
                   ipiv.data());
  } else if (is_complex_float<T>::value) {
    info = LAPACKE_cgetrf(LAPACK_ROW_MAJOR, m, n,
                          (lapack_complex_float *)b.data(), lda, ipiv.data());
    LAPACKE_cgetri(LAPACK_ROW_MAJOR, n, (lapack_complex_float *)b.data(), lda,
                   ipiv.data());
  } else {
    err_msg("inverse(): unspported element type.");
  }

  if (info) return false;

  return true;
}

template <typename T>
inline Matrix<T, 2> inverse(const Matrix<T, 2> &a) {
  assert(a.n_rows() == a.n_cols());

  Matrix<T, 2> res;
  inverse(res, a);

  return res;
}

template <typename T>
inline Matrix<T, 2> inverse(const MatrixRef<T, 2> &a) {
  assert(a.n_rows() == a.n_cols());

  Matrix<T, 2> res;
  inverse(res, a);

  return res;
}

template <typename T>
inline bool inv(Matrix<T, 2> &b, const Matrix<T, 2> &a) {
  return inverse(b, a);
}

template <typename T>
inline bool inv(Matrix<T, 2> &b, const MatrixRef<T, 2> &a) {
  return inverse(b, a);
}

template <typename T>
inline Matrix<T, 2> inv(const Matrix<T, 2> &a) {
  return inverse(a);
}

template <typename T>
inline Matrix<T, 2> inv(const MatrixRef<T, 2> &a) {
  return inverse(a);
}

}  // namespace slab

#endif  // SLAB_MATRIX_FNS_INV_H_
