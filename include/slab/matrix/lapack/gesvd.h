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

/// @file gesvd.h
/// @brief

#ifndef _SLAB_MATRIX_LAPACK_GESVD_H
#define _SLAB_MATRIX_LAPACK_GESVD_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup lapack_interface LAPACK INTERFACE
/// @{

template <typename T>
inline int lapack_gesvd(char jobu, char jobvt, Matrix<T, 2> &a, Matrix<T, 1> &s,
                        Matrix<T, 2> &u, Matrix<T, 2> &vt,
                        Matrix<T, 1> &superb) {
  int m = a.n_rows();
  int n = a.n_cols();
  int lda = n;
  int ldu = m;
  int ldvt = n;

  s = Matrix<T, 1>(std::min(m, n));
  u = Matrix<T, 2>(m, m);
  vt = Matrix<T, 2>(n, n);
  superb = Matrix<T, 1>(std::min(m, n));

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

/// @}

_SLAB_END_NAMESPACE

#endif
