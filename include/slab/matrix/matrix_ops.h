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
#include <type_traits>
#include "slab/matrix/matrix.h"

#include <iostream>
using namespace std;

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

template<typename T, std::size_t N>
Matrix<T, N> operator*(const Matrix<T, N> &a, const Matrix<T, N> &b) {
  Matrix<T, N> res = a;
  res *= b;
  return res;
}

template<typename T, std::size_t N>
Matrix<T, N> operator/(const Matrix<T, N> &a, const Matrix<T, N> &b) {
  Matrix<T, N> res = a;
  res /= b;
  return res;
}

template<typename T>
Matrix<T, 1> matmul(const Matrix<T, 2> &a, const Matrix<T, 1> &b) {
  assert(a.extent(1) == b.extent(0));
  const std::size_t nr = a.extent(0);
  Matrix<T, 1> res(nr);

  if (std::is_floating_point<T>::value)
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                a.rows(), a.cols(), 1.0, a.data(), b.size(), b.data(), 1, 0.0, res.data(), 1);

  return res;
}

template<typename T>
Matrix<T, 2> matmul(const Matrix<T, 2> &a, const Matrix<T, 2> &b) {
  assert(a.extent(1) == b.extent(0));
  const std::size_t nrows = a.rows();
  const std::size_t ncols = b.cols();
  Matrix<T, 2> res(nrows, ncols);

  if (std::is_floating_point<T>::value)
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                a.rows(), b.cols(), a.cols(),
                1.0, a.data(), a.cols(),
                b.data(), b.cols(),
                0.0, res.data(), res.cols());

  return res;
}

template<typename T, std::size_t N, typename... Args>
auto reshape(const Matrix<T, N> &a, Args... args) -> decltype(Matrix<T, sizeof...(args)>()) {
  Matrix<T, sizeof...(args)> res(args...);
  cblas_dcopy(a.size(), a.data(), 1, res.data(), 1);

  return res;
}

template<typename T>
Matrix<T, 2> transpose(const Matrix<T, 2> &a) {
  Matrix<T, 2> res(a.cols(), a.rows());
  for (std::size_t i = 0; i < a.rows(); ++i) {
    for (std::size_t j = 0; j < a.rows(); ++j) {
      res(j, i) = a(i, j);
    }
  }

  return res;
}

//template<typename T>
//Matrix<T, 2> operator*(const Matrix<T, 1> &u, const Matrix<T, 1> &v) {
//  const std::size_t n = u.extent(0);
//  const std::size_t m = v.extent(0);
//  Matrix<T, 2> res(n, m);  // an n-by-m matrix
//  for (std::size_t i = 0; i != n; ++i)
//    for (std::size_t j = 0; j != m; ++j)
//      res(i, j) = u[i] * v[j];
//
//  return res;
//}
//
//template<typename T>
//Matrix<T, 1> operator*(const Matrix<T, 2> &m, const Matrix<T, 1> &v) {
//  assert(m.extent(1) == v.extent(0));
//
//  const std::size_t nr = m.extent(0);
//  const std::size_t nc = m.extent(1);
//  Matrix<T, 1> res(nr);
//  for (std::size_t i = 0; i != nr; ++i)
//    for (std::size_t j = 0; j != nc; ++j)
//      res(i) += m(i, j) * v(j);
//
//  return res;
//}
//
//template<typename T>
//Matrix<T, 2> operator*(const Matrix<T, 2> &m1, const Matrix<T, 2> &m2) {
//  const std::size_t n = m1.extent(0);
//  const std::size_t m = m1.extent(1);
//  assert(m == m2.extent(0));
//  const std::size_t p = m2.extent(1);
//
//  Matrix<T, 2> res(n, p);
//  for (std::size_t i = 0; i != n; ++i)
//    for (std::size_t j = 0; j != p; ++j)
//      for (std::size_t k = 0; k != m; ++k)
//        res(i, j) += m1(i, k) * m2(k, j);
//
//  return res;
//}
//
//template<typename T, typename U = T>
//T dot_product(const MatrixRef<T, 1> &a, const MatrixRef<U, 1> &b) {
//  return std::inner_product(a.begin(), a.end(), b.begin(), T{});
//}

#endif // SLAB_MATRIX_OPERATIONS_H_