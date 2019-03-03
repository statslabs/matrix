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

/// @file prod.h
/// @brief Product of array elements.

#ifndef SLAB_MATRIX_FNS_PROD_H_
#define SLAB_MATRIX_FNS_PROD_H_

namespace slab {

template <typename T>
inline T prod(const Matrix<T, 1> &x) {
  return std::accumulate(x.begin(), x.end(), T{1}, std::multiplies<T>());
}

template <typename T>
inline T prod(const MatrixRef<T, 1> &x) {
  return std::accumulate(x.begin(), x.end(), T{1}, std::multiplies<T>());
}

template <typename T>
inline Matrix<T, 1> prod(const Matrix<T, 2> &x) {
  Matrix<T, 1> res(x.n_cols());
  for (std::size_t i = 0; i != x.n_cols(); ++i) {
    Matrix<T, 1> xcol = x.col(i);
    res(i) =
        std::accumulate(xcol.begin(), xcol.end(), T{1}, std::multiplies<T>());
  }

  return res;
}

template <typename T>
inline Matrix<T, 1> prod(const MatrixRef<T, 2> &x) {
  Matrix<T, 1> res(x.n_cols());
  for (std::size_t i = 0; i != x.n_cols(); ++i) {
    Matrix<T, 1> xcol = x.col(i);
    res(i) =
        std::accumulate(xcol.begin(), xcol.end(), T{1}, std::multiplies<T>());
  }

  return res;
}

}  // namespace slab

#endif  // SLAB_MATRIX_FNS_PROD_H_
