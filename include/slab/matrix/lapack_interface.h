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

#include <cassert>
#include <cstddef>

#include <complex>

#ifdef USE_MKL
#include "mkl.h"
#else
extern "C" {
#include "lapacke.h"
}
#endif

#include "slab/matrix/error.h"
#include "slab/matrix/matrix.h"
#include "slab/matrix/matrix_base.h"
#include "slab/matrix/traits.h"

#include "slab/matrix/lapack/gesv.h"

namespace slab {

/// @addtogroup lapack_interface LAPACK INTERFACE
/// @{

/// @addtogroup lapack_linear_equation_routines LAPACK Linear Equation Routines
/// @{

template <typename T>
inline int lapack_getrf(Matrix<T, 2> &a, Matrix<int, 1> &ipiv) {
  int info = 0;

  const int m = a.n_rows();
  const int n = a.n_cols();

  // assert(ipiv.size() >= std::max(1, std::min(m, n)));
  ipiv.clear();
  ipiv = Matrix<int, 1>(std::max(1, std::min(m, n)));

  const int lda = n;

  if (is_double<T>::value) {
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, (double *)a.data(), lda,
                          ipiv.data());
  } else if (is_float<T>::value) {
    info = LAPACKE_sgetrf(LAPACK_ROW_MAJOR, m, n, (float *)a.data(), lda,
                          ipiv.data());
  }

  return info;
}

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

template <typename T>
inline int lapack_gesvd(char jobu, char jobvt, Matrix<T, 2> &a, Matrix<T, 1> &s,
                        Matrix<T, 2> &u, Matrix<T, 2> &vt,
                        Matrix<T, 1> &superb) {
  int m = a.n_rows();
  int n = a.n_cols();
  int lda = n;
  int ldu = m;
  int ldvt = n;

  s = zeros<Matrix<T, 1>>(std::min(m, n));
  u = zeros<Matrix<T, 2>>(m, m);
  vt = zeros<Matrix<T, 2>>(n, n);
  superb = zeros<Matrix<T, 1>>(std::min(m, n));

  int info = 0;
  if (is_double<T>::value) {
    info =
        LAPACKE_dgesvd(LAPACK_ROW_MAJOR, jobu, jobvt, m, n, (double *)a.data(),
                       lda, (double *)s.data(), (double *)u.data(), ldu,
                       (double *)vt.data(), ldvt, (double *)superb.data());
  } else if (is_float<T>::value) {
    info =
        LAPACKE_sgesvd(LAPACK_ROW_MAJOR, jobu, jobvt, m, n, (float *)a.data(),
                       lda, (float *)s.data(), (float *)u.data(), ldu,
                       (float *)vt.data(), ldvt, (float *)superb.data());
  }

  if (jobu == 'S') u = u.cols(0, std::min(m, n) - 1);
  if (jobvt == 'S') vt = vt.rows(0, std::min(m, n) - 1);

  return info;
}

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

}  // namespace slab

#endif  // SLAB_MATRIX_LAPACK_INTERFACE_H_
