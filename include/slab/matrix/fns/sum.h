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

/// @file sum.h
/// @brief Sum of array elements.

#ifndef STATSLABS_MATRIX_SUM_H_
#define STATSLABS_MATRIX_SUM_H_

namespace slab {

template <typename T, std::size_t N>
inline T sum(const Matrix<T, N> &x) {
  return std::accumulate(x.begin(), x.end(), T{0});
}

template <typename T, std::size_t N>
inline T sum(const MatrixRef<T, N> &x) {
  return std::accumulate(x.begin(), x.end(), T{0});
}

template <typename T>
inline Matrix<T, 1> sum(const Matrix<T, 2> &x) {
  Matrix<T, 1> res(x.n_cols());
  for (std::size_t i = 0; i != x.n_cols(); ++i) {
    Matrix<T, 1> xcol = x.col(i);
    res(i) = std::accumulate(xcol.begin(), xcol.end(), T{0});
  }

  return res;
}

template <typename T>
inline Matrix<T, 1> sum(const MatrixRef<T, 2> &x) {
  Matrix<T, 1> res(x.n_cols());
  for (std::size_t i = 0; i != x.n_cols(); ++i) {
    Matrix<T, 1> xcol = x.col(i);
    res(i) = std::accumulate(xcol.begin(), xcol.end(), T{0});
  }

  return res;
}

}  // namespace slab

#endif  // STATSLABS_MATRIX_SUM_H_
