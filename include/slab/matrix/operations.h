//
// Copyright 2018 The StatsLabs Authors.
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
// -----------------------------------------------------------------------------
// operations.h
// -----------------------------------------------------------------------------
//
#ifndef SLAB_MATRIX_OPERATIONS_H_
#define SLAB_MATRIX_OPERATIONS_H_

#include <cstddef>
#include <algorithm>
#include "slab/matrix/matrix.h"


template<typename T, std::size_t N>
Matrix<T, N> operator+(const Matrix<T, N> &a, const Matrix<T, N> &b) {
  Matrix<T, N> res = a;
  res += b;
  return res;
}

template<typename T, std::size_t N>
Matrix<T, N> operator-(const Matrix<T, N> &a, const Matrix<T, N> &b) {
  Matrix<T, N> res = a;
  res -= b;
  return res;
}

template<typename T>
Matrix<T, 2> operator*(const Matrix<T, 1> &u, const Matrix<T, 1> &v) {
  const std::size_t n = u.extent(0);
  const std::size_t m = v.extent(0);
  Matrix<T, 2> res(n, m);  // an n-by-m matrix
  for (std::size_t i = 0; i != n; ++i)
    for (std::size_t j = 0; j != m; ++j)
      res(i, j) = u[i] * v[j];

  return res;
}

template<typename T>
Matrix<T, 1> operator*(const Matrix<T, 2> &m, const Matrix<T, 1> &v) {
  assert(m.extent(1) == v.extent(0));

  const std::size_t nr = m.extent(0);
  const std::size_t nc = m.extent(1);
  Matrix<T, 1> res(nr);
  for (std::size_t i = 0; i != nr; ++i)
    for (std::size_t j = 0; j != nc; ++j)
      res(i) += m(i, j) * v(j);

  return res;
}

template<typename T>
Matrix<T, 2> operator*(const Matrix<T, 2> &m1, const Matrix<T, 2> &m2) {
  const std::size_t n = m1.extent(0);
  const std::size_t m = m1.extent(1);
  assert(m == m2.extent(0));
  const std::size_t p = m2.extent(1);

  Matrix<T, 2> res(n, p);
  for (std::size_t i = 0; i != n; ++i)
    for (std::size_t j = 0; j != p; ++j)
      for (std::size_t k = 0; k != m; ++k)
        res(i, j) += m1(i, k) * m2(k, j);

  return res;
}

template<typename T, typename U = T>
T dot_product(const MatrixRef<T, 1> &a, const MatrixRef<U, 1> &b) {
  return std::inner_product(a.begin(), a.end(), b.begin(), T{});
}

#endif // SLAB_MATRIX_OPERATIONS_H_