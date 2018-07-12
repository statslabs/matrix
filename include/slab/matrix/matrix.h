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
// matrix.h
// -----------------------------------------------------------------------------
//
#ifndef SLAB_MATRIX_MATRIX_H_
#define SLAB_MATRIX_MATRIX_H_

#include "slab/matrix/traits.h"
#include "slab/matrix/matrix_ref.h"

template<typename T, std::size_t N>
class Matrix {
 public:
  static constexpr std::size_t order_ = N;     // number of dimensions
  using value_type = T;
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;

  Matrix() = default;
  Matrix(Matrix &&) = default;                 // move
  Matrix &operator=(Matrix &&) = default;
  Matrix(Matrix const &) = default;            // copy
  Matrix &operator=(Matrix const &) = default;
  ~Matrix() = default;

  template<typename U>
  Matrix(const MatrixRef<U, N> &);             // construct from MatrixRef
  template<typename U>
  Matrix &operator=(const MatrixRef<U, N> &);  // assign from MatrixRef

  template<typename... Exts>
  explicit Matrix(Exts... exts);               // specify the extents

  Matrix(MatrixInitializer<T, N>);             // initialize from list
  Matrix &operator=(MatrixInitializer<T, N>);  // assign from list

  template<typename U,
      std::size_t NN = N,
      typename = Enable_if<(NN > 1)>,
      typename = Enable_if<Convertible<U, std::size_t>()>>
  Matrix(std::initializer_list<U>) = delete;
  template<typename U,
      std::size_t NN = N,
      typename = Enable_if<(NN > 1)>,
      typename = Enable_if<Convertible<U, std::size_t>()>>
  Matrix &operator=(std::initializer_list<U>) = delete;

  // number of dimensions
  static constexpr std::size_t order() { return order_; }
  // #elements in the nth dimension
  std::size_t extent(std::size_t n) const { assert(n < order_); return desc_.extents[n]; }
  // total number of elements
  std::size_t size() const { return elems_.size(); }
  // the slice defining subscripting
  const MatrixSlice<N> &descriptor() const { return desc_; }

  T *data() { return elems_.data(); }          // "flat" element access
  const T *data() const { return elems_.data(); }

  std::size_t rows() const { return desc_.extents[0]; }
  std::size_t cols() const { return desc_.extents[1]; }

  // m(i,j,k) subscripting with integers
  template<typename... Args>
  Enable_if<matrix_impl::Requesting_element<Args...>(), T &>
  operator()(Args... args);

  template<typename... Args>
  Enable_if<matrix_impl::Requesting_element<Args...>(), const T &>
  operator()(Args... args) const;

  // m(s1, s2, s3) subscripting with slides
  template<typename... Args>
  Enable_if<matrix_impl::Requesting_slice<Args...>(), MatrixRef<T, N>>
  operator()(const Args &... args);

  template<typename... Args>
  Enable_if<matrix_impl::Requesting_slice<Args...>(), const MatrixRef<T, N>>
  operator()(const Args &... args) const;

  // m[i] row access
  MatrixRef<T, N - 1> operator[](std::size_t i) { return row(i); }
  MatrixRef<const T, N - 1> operator[](std::size_t i) const { return row(i); }

  // row access
  MatrixRef<T, N - 1> row(std::size_t n);
  MatrixRef<const T, N - 1> row(std::size_t n) const;

  // column access
  MatrixRef<T, N - 1> col(size_t n);
  MatrixRef<const T, N - 1> col(size_t n) const;

  template<typename F>
  Matrix &apply(F f);                          // f(x) for every element x

  // f(x, mx) for corresponding elements of *this and m
  template<typename M, typename F>
  Enable_if<Matrix_type<M>(), Matrix &>
  apply(const M &m, F f);

  Matrix &operator=(const T &value);           // assignment with scalar

  Matrix &operator+=(const T &value);          // scalar addition
  Matrix &operator-=(const T &value);          // scalar subtraction
  Matrix &operator*=(const T &value);          // scalar multiplication
  Matrix &operator/=(const T &value);          // scalar division
  Matrix &operator%=(const T &value);          // scalar modulo

  // matrix addition
  template<typename M>
  Enable_if<Matrix_type<M>(), Matrix &> operator+=(const M &x);
  // matrix subtraction
  template<typename M>
  Enable_if<Matrix_type<M>(), Matrix &> operator-=(const M &x);
  template<typename M>
  Enable_if<Matrix_type<M>(), Matrix &> operator*=(const M &x);
  template<typename M>
  Enable_if<Matrix_type<M>(), Matrix &> operator/=(const M &x);
  template<typename M>
  Enable_if<Matrix_type<M>(), Matrix &> operator%=(const M &x);

  iterator begin() { return elems_.begin(); }
  const_iterator begin() const { return elems_.cbegin(); }
  iterator end() { return elems_.end(); }
  const_iterator end() const { return elems_.cend(); }

  void clear();

 private:
  MatrixSlice<N> desc_;   // slice defining extents in the N dimensions
  std::vector<T> elems_;  // the elements
};

template<typename T, std::size_t N>
template<typename U>
Matrix<T, N>::Matrix(const MatrixRef<U, N> &x)  // copy desc_ and elements
    :desc_{x.descriptor().extents}, elems_{x.begin(), x.end()} {
  static_assert(Convertible<U, T>(),
                "Matrix constructor: incompatible element types");
}

template<typename T, std::size_t N>
template<typename U>
Matrix<T, N> &Matrix<T, N>::operator=(const MatrixRef<U, N> &x) {
  static_assert(Convertible<U, T>(), "Matrix =: incompatible element types");

  desc_ = x.descriptor();
  elems_.assign(x.begin(), x.end());
  return *this;
}

template<typename T, std::size_t N>
template<typename... Exts>
Matrix<T, N>::Matrix(Exts... exts)
    :desc_{exts...},     // copy extents
     elems_(desc_.size)  // allocate desc_.size elements and default initialize them
{}

template<typename T, std::size_t N>
Matrix<T, N>::Matrix(MatrixInitializer<T, N> init) {
  desc_.extents = matrix_impl::derive_extents<N>(init);
  desc_.size = matrix_impl::compute_strides(desc_.extents, desc_.strides);
  elems_.reserve(desc_.size);              // make room for slices
  matrix_impl::insert_flat(init, elems_);  // initialize from initializer list
  assert(elems_.size() == desc_.size);
}

template<typename T, std::size_t N>
template<typename... Args>
Enable_if<matrix_impl::Requesting_element<Args...>(), T &>
Matrix<T, N>::operator()(Args... args) {
  assert(matrix_impl::check_bounds(desc_, args...));
  return *(data() + desc_(args...));
}

template<typename T, std::size_t N>
template<typename... Args>
Enable_if<matrix_impl::Requesting_element<Args...>(), const T &>
Matrix<T, N>::operator()(Args... args) const {
  assert(matrix_impl::check_bounds(desc_, args...));
  return *(data() + desc_(args...));
}

template<typename T, std::size_t N>
template<typename... Args>
Enable_if<matrix_impl::Requesting_slice<Args...>(), MatrixRef<T, N>>
Matrix<T, N>::operator()(const Args &... args) {
  MatrixSlice<N> d;
  d.start = matrix_impl::do_slice(desc_, d, args...);
  return {d, data()};
}

template<typename T, std::size_t N>
template<typename... Args>
Enable_if<matrix_impl::Requesting_slice<Args...>(), const MatrixRef<T, N>>
Matrix<T, N>::operator()(const Args &... args) const {
  MatrixSlice<N> d;
  d.start = matrix_impl::do_slice(desc_, d, args...);
  return {d, data()};
}

// row
template<typename T, std::size_t N>
MatrixRef<T, N - 1> Matrix<T, N>::row(std::size_t n) {
  assert(n < rows());
  MatrixSlice<N - 1> row;
  matrix_impl::slice_dim<0>(n, desc_, row);
  return {row, data()};
}

template<typename T, std::size_t N>
MatrixRef<const T, N - 1> Matrix<T, N>::row(std::size_t n) const {
  assert(n < rows());
  MatrixSlice<N - 1> row;
  matrix_impl::slice_dim<0>(n, desc_, row);
  return {row, data()};
}

// col
template<typename T, std::size_t N>
MatrixRef<T, N - 1> Matrix<T, N>::col(std::size_t n) {
  assert(n < cols());
  MatrixSlice<N - 1> col;
  matrix_impl::slice_dim<1>(n, desc_, col);
  return {col, data()};
}

template<typename T, std::size_t N>
MatrixRef<const T, N - 1> Matrix<T, N>::col(std::size_t n) const {
  assert(n < cols());
  MatrixSlice<N - 1> col;
  matrix_impl::slice_dim<1>(n, desc_, col);
  return {col, data()};
}

template<typename T, std::size_t N>
template<typename F>
Matrix<T, N> &Matrix<T, N>::apply(F f) {
  for (auto &x : elems_) f(x);
  return *this;
}

template<typename T, std::size_t N>
template<typename M, typename F>
Enable_if<Matrix_type<M>(), Matrix<T, N> &> Matrix<T, N>::apply(const M &m, F f) {
  assert(same_extents(desc_, m.descriptor()));
  auto j = m.begin();
  for (auto i = begin(); i != end(); ++i) {
    f(*i, *j);
    ++j;
  }

  return *this;
}

template<typename T, std::size_t N>
Matrix<T, N> &Matrix<T, N>::operator+=(const T &val) {
  return apply([&](T &a) { a += val; });
}

template<typename T, std::size_t N>
Matrix<T, N> &Matrix<T, N>::operator-=(const T &val) {
  return apply([&](T &a) { a -= val; });
}

template<typename T, std::size_t N>
Matrix<T, N> &Matrix<T, N>::operator*=(const T &val) {
  return apply([&](T &a) { a *= val; });
}

template<typename T, std::size_t N>
Matrix<T, N> &Matrix<T, N>::operator/=(const T &val) {
  return apply([&](T &a) { a /= val; });
}

template<typename T, std::size_t N>
Matrix<T, N> &Matrix<T, N>::operator%=(const T &val) {
  return apply([&](T &a) { a %= val; });
}

template<typename T, std::size_t N>
template<typename M>
Enable_if<Matrix_type<M>(), Matrix<T, N> &> Matrix<T, N>::operator+=(const M &m) {
  //static_assert(m.order_ == N, "+=: mismatched Matrix dimensions");
  assert(same_extents(desc_, m.descriptor()));  // make sure sizes match

  return apply(m, [&](T &a, const Value_type<M> &b) { a += b; });
}

template<typename T, std::size_t N>
template<typename M>
Enable_if<Matrix_type<M>(), Matrix<T, N> &> Matrix<T, N>::operator-=(const M &m) {
  //static_assert(m.order_ == N, "+=: mismatched Matrix dimensions");
  assert(same_extents(desc_, m.descriptor()));  // make sure sizes match

  return apply(m, [&](T &a, const Value_type<M> &b) { a -= b; });
}

template<typename T, std::size_t N>
template<typename M>
Enable_if<Matrix_type<M>(), Matrix<T, N> &> Matrix<T, N>::operator*=(const M &m) {
  assert(same_extents(desc_, m.descriptor()));  // make sure sizes match

  return apply(m, [&](T &a, const Value_type<M> &b) { a *= b; });
}

template<typename T, std::size_t N>
template<typename M>
Enable_if<Matrix_type<M>(), Matrix<T, N> &> Matrix<T, N>::operator/=(const M &m) {
  assert(same_extents(desc_, m.descriptor()));  // make sure sizes match

  return apply(m, [&](T &a, const Value_type<M> &b) { a /= b; });
}

template<typename T, std::size_t N>
template<typename M>
Enable_if<Matrix_type<M>(), Matrix<T, N> &> Matrix<T, N>::operator%=(const M &m) {
  assert(same_extents(desc_, m.descriptor()));  // make sure sizes match

  return apply(m, [&](T &a, const Value_type<M> &b) { a %= b; });
}

template<typename T, std::size_t N>
void Matrix<T, N>::clear() {
  desc_.clear();
  elems_.clear();
}

template<typename T>
class Matrix<T, 0> {
 public:
  static constexpr std::size_t order_ = 0;
  using value_type = T;

  Matrix(const T &x = T{}) : elem_(x) {}

  Matrix &operator=(const T &value) {
    elem_ = value;
    return *this;
  }

  T &operator()() { return elem_; }

  const T &operator()() const { return elem_; }

  const MatrixSlice<0> &descriptor() const { return desc_; }

 private:
  MatrixSlice<0> desc_;
  T elem_;
};

////////////////////////////////////////
/// PRINTING UTILS
///////////////////////////////////////

// print Matrix, MatrixRef
template<typename T, std::size_t N>
std::ostream &operator<<(std::ostream &os, const Matrix<T, N> &m) {
  os << std::endl << '{';
  for (size_t i = 0; i != m.rows(); ++i) {
    os << m[i];
    if (i + 1 != m.rows()) os << ',';
  }
  return os << '}' << std::endl;
}

template<typename T, std::size_t N>
std::ostream &operator<<(std::ostream &os, const MatrixRef<T, N> &m) {
  os << std::endl << '{';
  for (size_t i = 0; i != m.rows(); ++i) {
    os << m[i];
    if (i + 1 != m.rows()) os << ',';
  }
  return os << '}' << std::endl;
}

//template<typename M>
//Enable_if<Matrix_type<M>(), std::ostream &>
//operator<<(std::ostream &os, const M &m) {
//    os << '{';
//    for (size_t i = 0; i != m.rows(); ++i) {
//        os << m[i];
//        if (i + 1 != m.rows()) os << ',';
//    }
//    return os << '}';
//}

template<typename T>
std::ostream &operator<<(std::ostream &os, const Matrix<T, 0> &m0) {
  return os << (const T &) m0;
}

#endif // SLAB_MATRIX_MATRIX_H_
