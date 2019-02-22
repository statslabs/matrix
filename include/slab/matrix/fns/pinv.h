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

/// @file pinv.h
/// @brief pseudo-inverse.

#ifndef SLAB_MATRIX_FNS_PINV_H_
#define SLAB_MATRIX_FNS_PINV_H_

namespace slab {

template <typename T>
inline bool pinv(Matrix<T, 2> &a_inv, const Matrix<T, 2> &a) {
  int m = a.n_rows();
  int n = a.n_cols();
  int k = std::min(m, n);

  Matrix<T, 2> a_copy = a, u, vt;
  Matrix<T, 1> s, superb;

  int info = lapack_gesvd('S', 'S', a_copy, s, u, vt, superb);
  if (info) return false;

  for (int i = 0; i < k; i++) {
    double ss;
    if (s[i] > 1.0e-9)
      ss = 1.0 / s[i];
    else
      ss = s[i];

    Matrix<T, 1> ui = u.col(i);
    blas_scal(ss, ui);
    u.col(i) = ui;
  }

  a_inv = zeros<Matrix<T, 2>>(n, m);
  T alpha = 1.0;
  T beta = 0.0;
  blas_gemm(CblasTrans, CblasTrans, alpha, vt, u, beta, a_inv);

  return true;
}

}  // namespace slab

#endif  // SLAB_MATRIX_FNS_PINV_H_
