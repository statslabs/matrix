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

/// @file kron.h
/// @brief Kronecker tensor product.

#ifndef SLAB_MATRIX_FNS_KRON_H_
#define SLAB_MATRIX_FNS_KRON_H_

namespace slab {

template <typename T>
inline Matrix<T, 2> kron(const Matrix<T, 2> &a, const Matrix<T, 2> &b) {
  const std::size_t a_rows = a.n_rows();
  const std::size_t a_cols = a.n_cols();
  const std::size_t b_rows = b.n_rows();
  const std::size_t b_cols = b.n_cols();

  Matrix<T, 2> res(a_rows * b_rows, a_cols * b_cols);
  for (std::size_t j = 0; j != a_cols; ++j) {
    for (std::size_t i = 0; i != a_rows; ++i) {
      res(slice{i * b_rows, b_rows}, slice{j * b_cols, b_cols}) = a(i, j) * b;
    }
  }

  return res;
}

}  // namespace slab

#endif  // SLAB_MATRIX_FNS_KRON_H_
