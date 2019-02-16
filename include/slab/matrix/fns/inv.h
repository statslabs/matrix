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

#ifndef STATSLABS_MATRIX_INV_H_
#define STATSLABS_MATRIX_INV_H_

namespace slab {

template <typename T>
inline bool inv(Matrix<T, 2> &b, const Matrix<T, 2> &a) {
  return true;
}

template <>
inline bool inv(Matrix<double, 2> &b, const Matrix<double, 2> &a) {
  b = eye<Matrix<double, 2>>(a.n_rows(), a.n_cols());

  int n = a.n_rows();
  int nrhs = b.n_cols();
  int lda = a.n_cols();
  int ldb = b.n_cols();

  Matrix<double, 2> a_copy = a;
  Matrix<int, 1> ipiv(n);

  int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, (double *)a_copy.data(),
                           lda, ipiv.data(), (double *)b.data(), ldb);
  if (info) return false;

  return true;
}

}  // namespace slab

#endif  // STATSLABS_MATRIX_INV_H_
