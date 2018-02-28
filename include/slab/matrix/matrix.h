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
  static constexpr std::size_t order = N;
  using value_type = T;
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;

  Matrix() = default;
  Matrix(Matrix &&) = default;                                   // move
  Matrix &operator=(Matrix &&) = default;
  Matrix(Matrix const &) = default;                              // copy
  Matrix &operator=(Matrix const &) = default;
  ~Matrix() = default;

  template<typename U>
  Matrix(const MatrixRef<U, N> &);                               // construct from MatrixRef
  template<typename U>
  Matrix &operator=(const MatrixRef<U, N> &);                    // assign from MatrixRef

  template<typename... Exts>
  explicit Matrix(Exts... exts);                                 // specify the extents

  Matrix(MatrixInitializer<T, N>);                               // initialize from list
  Matrix &operator=(MatrixInitializer<T, N>);                    // assign from list

  template<typename U, std::size_t NN = N, typename = Enable_if<(NN > 1)>,
      typename = Enable_if<Convertible<U, std::size_t>()>>
  Matrix(std::initializer_list<U>) = delete;
  template<typename U, std::size_t NN = N, typename = Enable_if<(NN > 1)>,
      typename = Enable_if<Convertible<U, std::size_t>()>>
  Matrix &operator=(std::initializer_list<U>) = delete;

  std::size_t extent(size_t n) const { return desc.extents[n]; } // #elements in the nth dimension
  std::size_t size() const { return elems.size(); }              // total number of elements
  std::size_t rows() const { return desc.extents[0]; }
  std::size_t cols() const { return desc.extents[1]; }
  const MatrixSlice<N> &descriptor() const { return desc; }

  T *data() { return elems.data(); }                             // "flat" element access
  const T *data() const { return elems.data(); }

  template<typename... Args>
  Enable_if<matrix_impl::Requesting_element<Args...>(), T &>
  operator()(Args... args);

  template<typename... Args>
  Enable_if<matrix_impl::Requesting_element<Args...>(), const T &>
  operator()(Args... args) const;                                  // subscripting with integers

  template<typename... Args>
  Enable_if<matrix_impl::Requesting_slice<Args...>(), MatrixRef<T, N>>
  operator()(const Args &... args);

  template<typename... Args>
  Enable_if<matrix_impl::Requesting_slice<Args...>(), const MatrixRef<T, N>>
  operator()(const Args &... args) const;                          // subscripting with slides

  MatrixRef<T, N - 1> operator[](std::size_t i) { return row(i); } // m[i] row access
  MatrixRef<const T, N - 1> operator[](std::size_t i) const { return row(i); }

  MatrixRef<T, N - 1> row(std::size_t n);                          // row access
  MatrixRef<const T, N - 1> row(std::size_t n) const;

  MatrixRef<T, N - 1> col(size_t n);                               // column access
  MatrixRef<const T, N - 1> col(size_t n) const;

  template<typename F>
  Matrix &apply(F f);                 // f(x) for every element x

  template<typename M, typename F>
  Enable_if<Matrix_type<M>(), Matrix &>
  apply(const M &m, F f);             // f(x, m) for corresponding elements

  Matrix &operator=(const T &value);  // assignment with scalar

  Matrix &operator+=(const T &value); // scalar addition
  Matrix &operator-=(const T &value); // scalar subtraction
  Matrix &operator*=(const T &value); // scalar multiplication
  Matrix &operator/=(const T &value); // scalar division
  Matrix &operator%=(const T &value); // scalar modulo

  template<typename M>
  Enable_if<Matrix_type<M>(), Matrix &> operator+=(const M &x);    // matrix addtion
  template<typename M>
  Enable_if<Matrix_type<M>(), Matrix &> operator-=(const M &x);    // matrix substraction

  iterator begin() { return elems.begin(); }
  const_iterator begin() const { return elems.cbegin(); }
  iterator end() { return elems.end(); }
  const_iterator end() const { return elems.cend(); }

 private:
  MatrixSlice<N> desc;
  std::vector<T> elems;
};

template<typename T, std::size_t N>
template<typename U>
Matrix<T, N>::Matrix(const MatrixRef<U, N> &x)
    :desc(x.descriptor().extents), elems(x.begin(), x.end()) // copy desc and elements
{
  static_assert(Convertible<U, T>(),
                "Matrix constructor: incompatible element types");
}

template<typename T, std::size_t N>
template<typename U>
Matrix<T, N> &Matrix<T, N>::operator=(const MatrixRef<U, N> &x) {
  static_assert(Convertible<U, T>(), "Matrix =: incompatible element types");

  desc = x.descriptor().extents;
  elems.assign(x.begin(), x.end());

  return *this;
}

template<typename T, std::size_t N>
template<typename... Exts>
Matrix<T, N>::Matrix(Exts... exts)
    : desc(exts...),    // copy extents
      elems(desc.size)  // allocate desc.size elements and default initialize them
{}

template<typename T, std::size_t N>
Matrix<T, N>::Matrix(MatrixInitializer<T, N> init) {
  desc.extents = matrix_impl::derive_extents<N>(init);
  desc.size = matrix_impl::compute_strides(desc.extents, desc.strides);
  elems.reserve(desc.size);
  matrix_impl::insert_flat(init, elems);
  assert(elems.size() == desc.size);
}

template<typename T, std::size_t N>
template<typename... Args>
Enable_if<matrix_impl::Requesting_element<Args...>(), T &>
Matrix<T, N>::operator()(Args... args) {
  assert(matrix_impl::check_bounds(desc, args...));
  return *(data() + desc(args...));
}

template<typename T, std::size_t N>
template<typename... Args>
Enable_if<matrix_impl::Requesting_element<Args...>(), const T &>
Matrix<T, N>::operator()(Args... args) const {
  assert(matrix_impl::check_bounds(desc, args...));
  return *(data() + desc(args...));
}

template<typename T, std::size_t N>
template<typename... Args>
Enable_if<matrix_impl::Requesting_slice<Args...>(), MatrixRef<T, N>>
Matrix<T, N>::operator()(const Args &... args) {
  MatrixSlice<N> d;
  d.start = matrix_impl::do_slice(desc, d, args...);
  return {d, data()};
}

template<typename T, std::size_t N>
template<typename... Args>
Enable_if<matrix_impl::Requesting_slice<Args...>(), const MatrixRef<T, N>>
Matrix<T, N>::operator()(const Args &... args) const {
  MatrixSlice<N> d;
  d.start = matrix_impl::do_slice(desc, d, args...);
  return {d, data()};
}

// row
template<typename T, std::size_t N>
MatrixRef<T, N - 1> Matrix<T, N>::row(std::size_t n) {
  assert(n < rows());
  MatrixSlice<N - 1> row;
  matrix_impl::slice_dim<0>(n, desc, row);
  return {row, data()};
}

template<typename T, std::size_t N>
MatrixRef<const T, N - 1> Matrix<T, N>::row(std::size_t n) const {
  assert(n < rows());
  MatrixSlice<N - 1> row;
  matrix_impl::slice_dim<0>(n, desc, row);
  return {row, data()};
}

// col
template<typename T, std::size_t N>
MatrixRef<T, N - 1> Matrix<T, N>::col(std::size_t n) {
  assert(n < cols());
  MatrixSlice<N - 1> col;
  matrix_impl::slice_dim<1>(n, desc, col);
  return {col, data()};
}

template<typename T, std::size_t N>
MatrixRef<const T, N - 1> Matrix<T, N>::col(std::size_t n) const {
  assert(n < cols());
  MatrixSlice<N - 1> col;
  matrix_impl::slice_dim<1>(n, desc, col);
  return {col, data()};
}

template<typename T, std::size_t N>
template<typename F>
Matrix<T, N> &Matrix<T, N>::apply(F f) {
  for (auto &x : elems) f(x);
  return *this;
}

template<typename T, std::size_t N>
template<typename M, typename F>
Enable_if<Matrix_type<M>(), Matrix<T, N> &> Matrix<T, N>::apply(const M &m, F f) {
  assert(same_extents(desc, m.descriptor()));
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
  //static_assert(m.order == N, "+=: mismatched Matrix dimensions");
  assert(same_extents(desc, m.descriptor()));  // make sure sizes match

  return apply(m, [&](T &a, const Value_type<M> &b) { a += b; });
}

template<typename T, std::size_t N>
template<typename M>
Enable_if<Matrix_type<M>(), Matrix<T, N> &> Matrix<T, N>::operator-=(const M &m) {
  //static_assert(m.order == N, "+=: mismatched Matrix dimensions");
  assert(same_extents(desc, m.descriptor()));  // make sure sizes match

  return apply(m, [&](T &a, const Value_type<M> &b) { a -= b; });
}

template<typename T>
class Matrix<T, 0> {
 public:
  static constexpr std::size_t order = 0;
  using value_type = T;

  Matrix(const T &x = T{}) : elem(x) {}

  Matrix &operator=(const T &value) {
    elem = value;
    return *this;
  }

  T &operator()() { return elem; }

  const T &operator()() const { return elem; }

  const MatrixSlice<0> &descriptor() const { return desc; }

 private:
  MatrixSlice<0> desc;
  T elem;
};

//template<typename T>
//T &Matrix<T, 1>::row(std::size_t i) { return &elems[i]; }

//template<typename T>
//T &Matrix<T, 0>::row(std::size_t i) = delete;

// print Matrix, MatrixRef
template<typename T, std::size_t N>
std::ostream &operator<<(std::ostream &os, const Matrix<T, N> &m) {
  os << '{';
  for (size_t i = 0; i != m.rows(); ++i) {
    os << m[i];
    if (i + 1 != m.rows()) os << ',';
  }
  return os << '}';
}

template<typename T, std::size_t N>
std::ostream &operator<<(std::ostream &os, const MatrixRef<T, N> &m) {
  os << '{';
  for (size_t i = 0; i != m.rows(); ++i) {
    os << m[i];
    if (i + 1 != m.rows()) os << ',';
  }
  return os << '}';
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
