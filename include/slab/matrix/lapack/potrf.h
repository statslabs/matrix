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

/// @file potrf.h
/// @brief

#ifndef _SLAB_MATRIX_LAPACK_POTRF_H
#define _SLAB_MATRIX_LAPACK_POTRF_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup lapack_interface LAPACK INTERFACE
/// @{

template <typename T>
inline int lapack_potrf(Matrix<T, 2> &a) {
  char uplo = 'U';
  int n = a.n_rows();
  int lda = a.n_cols();

  int info = 0;
  if (is_double<T>::value) {
    info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, uplo, n, (double *)a.data(), lda);
  } else if (is_float<T>::value) {
    info = LAPACKE_spotrf(LAPACK_ROW_MAJOR, uplo, n, (float *)a.data(), lda);
  } else if (is_complex_double<T>::value) {
    info = LAPACKE_zpotrf(LAPACK_ROW_MAJOR, uplo, n,
                          reinterpret_cast<lapack_complex_double *>(a.data()),
                          lda);
  } else if (is_complex_float<T>::value) {
    info =
        LAPACKE_cpotrf(LAPACK_ROW_MAJOR, uplo, n,
                       reinterpret_cast<lapack_complex_float *>(a.data()), lda);
  }

  return info;
}

/// @}

_SLAB_END_NAMESPACE

#endif
