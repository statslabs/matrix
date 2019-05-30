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

/// @file gesv.h
/// @brief Computes the solution to the system of linear equations with a square
/// coefficient matrix A and multiple right-hand sides.

#ifndef _SLAB_MATRIX_LAPACK_GESV_H
#define _SLAB_MATRIX_LAPACK_GESV_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup lapack_interface LAPACK INTERFACE
/// @{

template <typename T>
inline int lapack_gesv(Matrix<T, 2> &a, Matrix<int, 1> &ipiv, Matrix<T, 2> &b) {
  assert(a.n_rows() == b.n_rows());

  int n = a.n_rows();
  int nrhs = b.n_cols();
  int lda = a.n_cols();
  int ldb = b.n_cols();

  int info = 0;
  if (is_double<T>::value) {
    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, (double *)a.data(), lda,
                         ipiv.data(), (double *)b.data(), ldb);
  } else if (is_float<T>::value) {
    info = LAPACKE_sgesv(LAPACK_ROW_MAJOR, n, nrhs, (float *)a.data(), lda,
                         ipiv.data(), (float *)b.data(), ldb);
  } else if (is_complex_double<T>::value) {
    info = LAPACKE_zgesv(
        LAPACK_ROW_MAJOR, n, nrhs,
        reinterpret_cast<lapack_complex_double *>(a.data()), lda, ipiv.data(),
        reinterpret_cast<lapack_complex_double *>(b.data()), ldb);
  } else if (is_complex_float<T>::value) {
    info = LAPACKE_cgesv(
        LAPACK_ROW_MAJOR, n, nrhs,
        reinterpret_cast<lapack_complex_float *>(a.data()), lda, ipiv.data(),
        reinterpret_cast<lapack_complex_float *>(b.data()), ldb);
  }

  return info;
}

/// @}

_SLAB_END_NAMESPACE

#endif
