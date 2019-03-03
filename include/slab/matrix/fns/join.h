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

/// @file join.h
/// @brief concatenation of matrices.

#ifndef SLAB_MATRIX_FNS_JOIN_H_
#define SLAB_MATRIX_FNS_JOIN_H_

namespace slab {

// join_vecs()
template <typename T>
inline Matrix<T, 1> join_vecs(const Matrix<T, 1> &a, const Matrix<T, 1> &b) {
  Matrix<T, 1> res(a.n_rows() + b.n_rows());
  res(slice{0, a.n_rows()}) = a;
  res(slice{a.n_rows(), b.n_rows()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 1> join_vecs(const MatrixRef<T, 1> &a,
                              const MatrixRef<T, 1> &b) {
  Matrix<T, 1> res(a.n_rows() + b.n_rows());
  res(slice{0, a.n_rows()}) = a;
  res(slice{a.n_rows(), b.n_rows()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 1> join_vecs(const Matrix<T, 1> &a, const MatrixRef<T, 1> &b) {
  Matrix<T, 1> res(a.n_rows() + b.n_rows());
  res(slice{0, a.n_rows()}) = a;
  res(slice{a.n_rows(), b.n_rows()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 1> join_vecs(const MatrixRef<T, 1> &a, const Matrix<T, 1> &b) {
  Matrix<T, 1> res(a.n_rows() + b.n_rows());
  res(slice{0, a.n_rows()}) = a;
  res(slice{a.n_rows(), b.n_rows()}) = b;

  return res;
}

template <typename T, typename... Args>
inline Matrix<T, 1> join_vecs(const Matrix<T, 1> &x, Args... args) {
  return join_vecs(x, join_vecs(args...));
}

template <typename T, typename... Args>
inline Matrix<T, 1> join_vecs(const MatrixRef<T, 1> &x, Args... args) {
  return join_vecs(x, join_vecs(args...));
}

// join_rows()

template <typename T>
inline Matrix<T, 2> join_rows(const Matrix<T, 2> &a, const Matrix<T, 2> &b) {
  Matrix<T, 2> res;
  if (a.empty() && b.empty())
    return res;
  else if (a.empty())
    return b;
  else if (b.empty())
    return a;
  else if (a.n_rows() != b.n_rows())
    err_quit("joint_rows(): inconsistent number of rows");

  res = slab::zeros<Matrix<T, 2>>(a.n_rows(), a.n_cols() + b.n_cols());
  res(slice{0, a.n_rows()}, slice{0, a.n_cols()}) = a;
  res(slice{0, a.n_rows()}, slice{a.n_cols(), b.n_cols()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 2> join_rows(const MatrixRef<T, 2> &a,
                              const MatrixRef<T, 2> &b) {
  assert(a.n_rows() == b.n_rows());

  Matrix<T, 2> res(a.n_rows(), a.n_cols() + b.n_cols());
  res(slice{0, a.n_rows()}, slice{0, a.n_cols()}) = a;
  res(slice{0, a.n_rows()}, slice{a.n_cols(), b.n_cols()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 2> join_rows(const Matrix<T, 2> &a, const MatrixRef<T, 2> &b) {
  assert(a.n_rows() == b.n_rows());

  Matrix<T, 2> res(a.n_rows(), a.n_cols() + b.n_cols());
  res(slice{0, a.n_rows()}, slice{0, a.n_cols()}) = a;
  res(slice{0, a.n_rows()}, slice{a.n_cols(), b.n_cols()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 2> join_rows(const MatrixRef<T, 2> &a, const Matrix<T, 2> &b) {
  assert(a.n_rows() == b.n_rows());

  Matrix<T, 2> res(a.n_rows(), a.n_cols() + b.n_cols());
  res(slice{0, a.n_rows()}, slice{0, a.n_cols()}) = a;
  res(slice{0, a.n_rows()}, slice{a.n_cols(), b.n_cols()}) = b;

  return res;
}

// join_cols()

template <typename T>
inline Matrix<T, 2> join_cols(const Matrix<T, 2> &a, const Matrix<T, 2> &b) {
  Matrix<T, 2> res;
  if (a.empty() && b.empty())
    return res;
  else if (a.empty())
    return b;
  else if (b.empty())
    return a;
  else if (a.n_cols() != b.n_cols())
    err_quit("joint_rows(): inconsistent number of columns");

  res = slab::zeros<Matrix<T, 2>>(a.n_rows() + b.n_rows(), a.n_cols());
  res(slice{0, a.n_rows()}, slice{0, a.n_cols()}) = a;
  res(slice{a.n_rows(), b.n_rows()}, slice{0, a.n_cols()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 2> join_cols(const MatrixRef<T, 2> &a,
                              const MatrixRef<T, 2> &b) {
  assert(a.n_cols() == b.n_cols());

  Matrix<T, 2> res(a.n_rows() + b.n_rows(), a.n_cols());
  res(slice{0, a.n_rows()}, slice{0, a.n_cols()}) = a;
  res(slice{a.n_rows(), b.n_rows()}, slice{0, a.n_cols()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 2> join_cols(const Matrix<T, 2> &a, const MatrixRef<T, 2> &b) {
  assert(a.n_cols() == b.n_cols());

  Matrix<T, 2> res(a.n_rows() + b.n_rows(), a.n_cols());
  res(slice{0, a.n_rows()}, slice{0, a.n_cols()}) = a;
  res(slice{a.n_rows(), b.n_rows()}, slice{0, a.n_cols()}) = b;

  return res;
}

template <typename T>
inline Matrix<T, 2> join_cols(const MatrixRef<T, 2> &a, const Matrix<T, 2> &b) {
  assert(a.n_cols() == b.n_cols());

  Matrix<T, 2> res(a.n_rows() + b.n_rows(), a.n_cols());
  res(slice{0, a.n_rows()}, slice{0, a.n_cols()}) = a;
  res(slice{a.n_rows(), b.n_rows()}, slice{0, a.n_cols()}) = b;

  return res;
}

}  // namespace slab

#endif  // SLAB_MATRIX_FNS_JOIN_H_
