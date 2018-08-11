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

/// @file triangular_matrix.h
/// @brief Triangular matrix

#ifndef STATSLABS_MATRIX_TRIANGULAR_MATRIX_H_
#define STATSLABS_MATRIX_TRIANGULAR_MATRIX_H_

#include "slab/matrix/packed_matrix.h"

template<typename T, typename TRI>
class TriangularMatrix : public PackedMatrix<T, TRI> {
 public:
  TriangularMatrix(std::size_t n) : PackedMatrix<T, TRI>{n} {}
  TriangularMatrix(std::size_t n, const Matrix<T, 1> &v) : PackedMatrix<T, TRI>(n, v) {}

  T &operator()(std::size_t i, std::size_t j) {
    if (!this->desc_.other_half(i, j))
      return *(this->data() + this->desc_(i, j));
    else
      return elem_in_other_half(i, j);
  }

  const T &operator()(std::size_t i, std::size_t j) const {
    if (!this->desc_.other_half(i, j))
      return *(this->data() + this->desc_(i, j));
    else
      return elem_in_other_half(i, j);
  }

  T &elem_in_other_half(std::size_t i, std::size_t j) {
    return T{0};
  }

  const T &elem_in_other_half(std::size_t i, std::size_t j) const {
    return T{0};
  }
};

template<typename T, typename TRI>
std::ostream &operator<<(std::ostream &os, const TriangularMatrix<T, TRI> &m) {
  for (std::size_t i = 0; i != m.n_rows(); ++i) {
    for (std::size_t j = 0; j != m.n_cols(); ++j) {
      os << m(i, j) << "\t";
    }
    os << std::endl;
  }

  return os << std::endl;
}

#endif // STATSLABS_MATRIX_TRIANGULAR_MATRIX_H_
