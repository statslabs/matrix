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

/// @file lu.h
/// @brief LU decomposition.

#ifndef _SLAB_MATRIX_FNS_TRANS_H
#define _SLAB_MATRIX_FNS_TRANS_H

_SLAB_BEGIN_NAMESPACE

template <typename T>
inline void lu(Matrix<T, 2> &L, Matrix<T, 2> &U, const Matrix<T, 2> &X)
{
  _SLAB_ASSERT(X.n_rows() == X.n_cols(), "X should be a square matrix");

  std::size_t n = X.n_rows();
  L = zeros<Matrix<T, 2>>(n, n);
  U = zeros<Matrix<T, 2>>(n, n);

  L.diag() = T{1};

  // Decompose matrix X into:
  // L -- Lower triangular matrix
  // U -- Upper triangular matrix
  for (std::size_t i = 0; i != n; ++i) {
    // Upper Triangular
    for (std::size_t k = i; k != n; ++k) {
      // Summation of L(i, j) * U(j, k)
      T sum = {};
      for (std::size_t j = 0; j < i; ++j) {
  	sum += L(i, j) * U(j, k);
      }

      // Evaluating U(i, k)
      U(i, k) = X(i, k) - sum;
    }

    // Lower Triangular
    for (std::size_t k = i + 1; k != n; ++k) {
      // Summation of L(k, j) * U(j, i)
      T sum = {};
      for (std::size_t j = 0; j < i; ++j) {
  	sum += L(k, j) * U(j, i);
      }

      // Evaluating L(k, i)
      L(k, i) = (X(k, i) - sum) / U(i, i);
    }
  }
}

_SLAB_END_NAMESPACE

#endif  // _SLAB_MATRIX_FNS_TRANS_H
