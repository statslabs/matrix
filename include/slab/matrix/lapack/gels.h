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

/// @file gels.h
/// @brief

#ifndef _SLAB_MATRIX_LAPACK_GELS_H
#define _SLAB_MATRIX_LAPACK_GELS_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup lapack_interface LAPACK INTERFACE
/// @{

template <typename T>
inline int lapack_gels(char trans, Matrix<T, 2> &a, Matrix<T, 2> &b) {
  int m = a.n_rows();
  int n = a.n_cols();
  int nrhs = b.n_cols();
  int lda = a.n_cols();
  int ldb = b.n_cols();

  int info = 0;
  if (is_double<T>::value) {
    info = LAPACKE_dgels(LAPACK_ROW_MAJOR, trans, m, n, nrhs,
                         (double *)a.data(), lda, (double *)b.data(), ldb);
  } else if (is_float<T>::value) {
    info = LAPACKE_sgels(LAPACK_ROW_MAJOR, trans, m, n, nrhs, (float *)a.data(),
                         lda, (float *)b.data(), ldb);
  } else if (is_complex_double<T>::value) {
    info =
        LAPACKE_zgels(LAPACK_ROW_MAJOR, trans, m, n, nrhs,
                      reinterpret_cast<lapack_complex_double *>(a.data()), lda,
                      reinterpret_cast<lapack_complex_double *>(b.data()), ldb);
  } else if (is_complex_float<T>::value) {
    info =
        LAPACKE_cgels(LAPACK_ROW_MAJOR, trans, m, n, nrhs,
                      reinterpret_cast<lapack_complex_float *>(a.data()), lda,
                      reinterpret_cast<lapack_complex_float *>(b.data()), ldb);
  }

  return info;
}

/// @}

_SLAB_END_NAMESPACE

#endif
