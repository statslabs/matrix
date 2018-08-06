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

/// @file lapack_interface.h
/// @brief LAPACK interface

#ifndef SLAB_MATRIX_LAPACK_INTERFACE_H_
#define SLAB_MATRIX_LAPACK_INTERFACE_H_

#include "slab/matrix/matrix.h"
#include "slab/matrix/traits.h"

/// @addtogroup lapack_interface LAPACK INTERFACE
/// @{

/// @addtogroup lapack_linear_equation_routines LAPACK Linear Equation Routines
/// @{

template<typename T>
int lapack_getrf(Matrix<T, 2> &a, Matrix<int, 1> &ipiv) {

  int info = 0;

  const int m = a.n_rows();
  const int n = a.n_cols();

  //assert(ipiv.size() >= std::max(1, std::min(m, n)));
  ipiv.clear();
  ipiv = Matrix<int, 1>(std::max(1, std::min(m, n)));

  const int lda = n;

  if (is_double<T>::value) {
    info = LAPACKE_dgetrf(
        LAPACK_ROW_MAJOR,
        m,
        n,
        (double *) a.data(),
        lda,
        ipiv.data()
    );
  } else if (is_float<T>::value) {
    info = LAPACKE_sgetrf(
        LAPACK_ROW_MAJOR,
        m,
        n,
        (float *) a.data(),
        lda,
        ipiv.data()
    );
  }

  return info;
}

template<typename T>
int lapack_gesv(Matrix<T, 2> &a, Matrix<int, 1> &ipiv, Matrix<T, 2> &b) {
  assert(a.n_rows() == b.n_rows());

  int n = a.n_rows();
  int nrhs = b.n_cols();
  int lda = a.n_cols();
  int ldb = b.n_cols();

  int info = 0;
  if (is_double<T>::value) {
    info = LAPACKE_dgesv(
        LAPACK_ROW_MAJOR,
        n,
        nrhs,
        (double *) a.data(),
        lda,
        ipiv.data(),
        (double *) b.data(),
        ldb
    );
  } else if (is_float<T>::value) {
    info = LAPACKE_sgesv(
        LAPACK_ROW_MAJOR,
        n,
        nrhs,
        (float *) a.data(),
        lda,
        ipiv.data(),
        (float *) b.data(),
        ldb
    );
  } else if (is_complex_double<T>::value) {
    info = LAPACKE_zgesv(
        LAPACK_ROW_MAJOR,
        n,
        nrhs,
        reinterpret_cast<lapack_complex_double *>(a.data()),
        lda,
        ipiv.data(),
        reinterpret_cast<lapack_complex_double *>(b.data()),
        ldb
    );
  } else if (is_complex_float<T>::value) {
    info = LAPACKE_cgesv(
        LAPACK_ROW_MAJOR,
        n,
        nrhs,
        reinterpret_cast<lapack_complex_float *>(a.data()),
        lda,
        ipiv.data(),
        reinterpret_cast<lapack_complex_float *>(b.data()),
        ldb
    );
  }

  return info;
}

/// @}
/// @}

#endif // SLAB_MATRIX_LAPACK_INTERFACE_H_
