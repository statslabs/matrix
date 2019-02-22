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

/// @file trans.h
/// @brief transpose of matrix.

#ifndef SLAB_MATRIX_FNS_TRANS_H_
#define SLAB_MATRIX_FNS_TRANS_H_

namespace slab {

template <typename T>
inline Matrix<T, 2> transpose(const Matrix<T, 1> &a) {
  Matrix<T, 2> res(1, a.n_rows());
  std::copy(a.cbegin(), a.cend(), res.begin());

  return res;
}

template <typename T>
inline Matrix<T, 2> transpose(const MatrixRef<T, 1> &a) {
  Matrix<T, 2> res(1, a.n_rows());
  std::copy(a.cbegin(), a.cend(), res.begin());

  return res;
}

template <typename T>
inline Matrix<T, 2> transpose(const Matrix<T, 2> &a) {
  Matrix<T, 2> res(a.n_cols(), a.n_rows());
  for (std::size_t i = 0; i < a.n_rows(); ++i) {
    res.col(i) = a.row(i);
  }

  return res;
}

template <typename T>
inline Matrix<T, 2> transpose(const MatrixRef<T, 2> &a) {
  Matrix<T, 2> res(a.n_cols(), a.n_rows());
  for (std::size_t i = 0; i < a.n_rows(); ++i) {
    res.col(i) = a.row(i);
  }

  return res;
}

}  // namespace slab

#endif  // SLAB_MATRIX_FNS_TRANS_H_
