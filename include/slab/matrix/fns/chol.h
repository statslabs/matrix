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

/// @file chol.h
/// @brief Cholesky decomposition.

#ifndef SLAB_MATRIX_FNS_CHOL_H_
#define SLAB_MATRIX_FNS_CHOL_H_

namespace slab {

template <typename T>
inline Matrix<T, 2> chol(const Matrix<T, 2> &x) {
  assert(x.n_rows() == x.n_cols());

  Matrix<T, 2> x_copy(x);
  int info = lapack_potrf(x_copy);
  if (info) _SLAB_ERROR("chol(): unsuccessful");

  Matrix<T, 2> res = zeros<Matrix<T, 2>>(x.n_rows(), x.n_cols());
  for (std::size_t i = 0; i != x.n_rows(); ++i) {
    for (std::size_t j = i; j != x.n_cols(); ++j) {
      res(i, j) = x_copy(i, j);
    }
  }

  return res;
}

}  // namespace slab

#endif  // SLAB_MATRIX_FNS_CHOL_H_
