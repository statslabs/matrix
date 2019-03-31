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

/// @file matrix_ref.h
/// @brief A MatrixRef template

#ifndef SLAB_MATRIX_MATRIX_REF_H_
#define SLAB_MATRIX_MATRIX_REF_H_

#include <cassert>
#include <cstddef>

#include <array>
#include <iterator>
#include <string>
#include <type_traits>

#include "slab/matrix/matrix_base.h"
#include "slab/matrix/matrix_slice.h"
#include "slab/matrix/packed_matrix.h"
#include "slab/matrix/support.h"

namespace slab {

template <typename T, std::size_t N>
class MatrixRefIterator;

template <typename T, std::size_t N>
class MatrixRef : public MatrixBase<T, N> {
  // ----------------------------------------
  // The core member functions in book 'TCPL'
  // ----------------------------------------
 public:
  //! @cond Doxygen_Suppress
  using iterator = MatrixRefIterator<T, N>;
  using const_iterator = MatrixRefIterator<const T, N>;

  MatrixRef() = delete;
  MatrixRef(MatrixRef &&) = default;  // move
  MatrixRef &operator=(MatrixRef &&);
  MatrixRef(const MatrixRef &) = default;  // copy
  MatrixRef &operator=(const MatrixRef &);
  ~MatrixRef() = default;
  //! @endcond

  //! construct from MatrixRef
  template <typename U>
  MatrixRef(const MatrixRef<U, N> &x);
  //! assign from MatrixRef
  template <typename U>
  MatrixRef &operator=(const MatrixRef<U, N> &x);

  //! construct from Matrix
  template <typename U>
  MatrixRef(const Matrix<U, N> &);
  //! assign from Matrix
  template <typename U>
  MatrixRef &operator=(const Matrix<U, N> &);

  //! assign from list
  MatrixRef &operator=(MatrixInitializer<T, N>);

  MatrixRef(const MatrixSlice<N> &s, T *p) : MatrixBase<T, N>{s}, ptr_{p} {}

  //! total number of elements
  std::size_t size() const { return this->desc_.size; }

  //! "flat" element access
  ///@{
  T *data() { return ptr_; }
  const T *data() const { return ptr_; }
  ///@}

 private:
  T *ptr_;

  // ---------------------------------------------
  // Member functions for subscripting and slicing
  // ---------------------------------------------
 public:
  //! m(i,j,k) subscripting with integers
  ///@{
  template <typename... Args>
  Enable_if<matrix_impl::Requesting_element<Args...>(), T &> operator()(
      Args... args) {
    return MatrixBase<T, N>::template operator()<Args...>(args...);
  }

  template <typename... Args>
  Enable_if<matrix_impl::Requesting_element<Args...>(), const T &> operator()(
      Args... args) const {
    return MatrixBase<T, N>::template operator()<Args...>(args...);
  }
  ///@}

  //! m(s1, s2, s3) subscripting with slides
  ///@{
  template <typename... Args>
  Enable_if<matrix_impl::Requesting_slice<Args...>(), MatrixRef<T, N>>
  operator()(const Args &... args);

  template <typename... Args>
  Enable_if<matrix_impl::Requesting_slice<Args...>(), MatrixRef<const T, N>>
  operator()(const Args &... args) const;
  ///@}

  //! m[i] row access
  ///@{
  MatrixRef<T, N - 1> operator[](std::size_t i) { return row(i); }
  MatrixRef<const T, N - 1> operator[](std::size_t i) const { return row(i); }
  ///@}

  //! row access
  ///@{
  MatrixRef<T, N - 1> row(std::size_t n);
  MatrixRef<const T, N - 1> row(std::size_t n) const;
  ///@}

  //! column access
  ///@{
  MatrixRef<T, N - 1> col(std::size_t n);
  MatrixRef<const T, N - 1> col(std::size_t n) const;
  ///@}

  //! multiple rows access
  ///@{
  MatrixRef<T, N> rows(std::size_t i, std::size_t j);
  MatrixRef<const T, N> rows(std::size_t i, std::size_t j) const;
  ///@}

  //! multiple columns access
  ///@{
  MatrixRef<T, N> cols(std::size_t i, std::size_t j);
  MatrixRef<const T, N> cols(std::size_t i, std::size_t j) const;
  ///@}

  //! element iterators
  ///@{
  iterator begin() { return {this->desc_, ptr_}; }
  const_iterator begin() const { return {this->desc_, ptr_}; }
  iterator end() { return {this->desc_, ptr_, true}; }
  const_iterator end() const { return {this->desc_, ptr_, true}; }
  ///@}

  // --------------------------------------------------
  // Member functions for matrix arithmetic operations
  // --------------------------------------------------
 public:
  //! matrix arithmetic operations
  ///@{
  template <typename F>
  MatrixRef &apply(F f);  // f(x) for every element x

  // f(x, mx) for corresponding elements of *this and m
  template <typename M, typename F>
  Enable_if<Matrix_type<M>(), MatrixRef &> apply(const M &m, F f);

  MatrixRef &operator=(const T &value);   // assignment with scalar
  MatrixRef &operator+=(const T &value);  // scalar addition
  MatrixRef &operator-=(const T &value);  // scalar subtraction
  MatrixRef &operator*=(const T &value);  // scalar multiplication
  MatrixRef &operator/=(const T &value);  // scalar division
  MatrixRef &operator%=(const T &value);  // scalar modulo

  // matrix addition
  template <typename M>
  Enable_if<Matrix_type<M>(), MatrixRef &> operator+=(const M &x);
  // matrix subtraction
  template <typename M>
  Enable_if<Matrix_type<M>(), MatrixRef &> operator-=(const M &x);
  // element-wise multiplication
  template <typename M>
  Enable_if<Matrix_type<M>(), MatrixRef &> operator*=(const M &x);
  // element-wise division
  template <typename M>
  Enable_if<Matrix_type<M>(), MatrixRef &> operator/=(const M &x);
  // element-wise modulus
  template <typename M>
  Enable_if<Matrix_type<M>(), MatrixRef &> operator%=(const M &x);

  template <typename U = typename std::remove_const<T>::type>
  Matrix<U, N> operator-() const;
  ///@}

  // -----------------------------------
  // Member functions for a MatrixRef<T, 1>
  // -----------------------------------

 public:
  //! Assign a vector from a matrix
  ///@{
  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 1)>>
    MatrixRef &operator=(const MatrixRef<U, 2> &x);
  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 1)>>
    MatrixRef &operator=(const Matrix<U, 2> &x);
  ///@}

  //! sub-vector access for Matrix<T, 1>
  ///@{
  template <std::size_t NN = N, typename = Enable_if<(NN == 1)>>
  MatrixRef<T, 1> subvec(std::size_t first_index, std::size_t last_index) {
    return this->operator()(slice{first_index, last_index - first_index + 1});
  }
  template <std::size_t NN = N, typename = Enable_if<(NN == 1)>>
  MatrixRef<const T, 1> subvec(std::size_t first_index,
                               std::size_t last_index) const {
    return this->operator()(slice{first_index, last_index - first_index + 1});
  }
  ///@}

  //! Construct a std::vector from a slab::MatrixRef<T, 1>
  template <std::size_t NN = N, typename = Enable_if<(NN == 1)>>
    std::vector<T> std_vec() const {
    return std::vector<T>(begin(), end());
  }

  // -----------------------------------
  // Member functions for a Matrix<T, 2>
  // -----------------------------------

 public:
  //! Assign a matrix from a vector
  ///@{
  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 2)>>
    MatrixRef &operator=(const MatrixRef<U, 1> &x);
  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 2)>>
    MatrixRef &operator=(const Matrix<U, 1> &x);
  ///@}

  //! sub-matrix access for Matrix<T, 2>
  ///@{
  template <std::size_t NN = N, typename = Enable_if<(NN == 2)>>
  MatrixRef<T, 2> submat(std::size_t first_row, std::size_t first_col,
                         std::size_t last_row, std::size_t last_col) {
    return this->operator()(slice{first_row, last_row - first_row + 1},
                            slice{first_col, last_col - first_col + 1});
  }

  template <std::size_t NN = N, typename = Enable_if<(NN == 2)>>
  MatrixRef<const T, 2> submat(std::size_t first_row, std::size_t first_col,
                               std::size_t last_row,
                               std::size_t last_col) const {
    return this->operator()(slice{first_row, last_row - first_row + 1},
                            slice{first_col, last_col - first_col + 1});
  }
  ///@}

  //! diagonal elements access for Matrix<T, 2>
  ///@{
  template <std::size_t NN = N, typename = Enable_if<(NN == 2)>>
  MatrixRef<T, 1> diag() {
    assert(this->n_rows() == this->n_cols());

    MatrixSlice<1> d;
    d.start = this->desc_.start;
    d.extents[0] = this->n_rows();
    d.strides[0] = this->n_rows() + 1;

    return {d, data()};
  }
  template <std::size_t NN = N, typename = Enable_if<(NN == 2)>>
  MatrixRef<const T, 1> diag() const {
    assert(this->n_rows() == this->n_cols());

    MatrixSlice<1> d;
    d.start = this->desc_.start;
    d.extents[0] = this->n_rows();
    d.strides[0] = this->n_rows() + 1;

    return {d, data()};
  }
  ///@}

  template <std::size_t NN = N, typename = Enable_if<(NN == 1) || (NN == 2)>>
  void print(const std::string &str = "") const {
    printf("\n %s\n", str.c_str());
    for (std::size_t i = 0; i != this->n_rows(); ++i) {
      raw_print(row(i));
      printf("\n");
    }
  }
};

template <typename T, std::size_t N>
MatrixRef<T, N> &MatrixRef<T, N>::operator=(MatrixRef &&x) {
  assert(same_extents(this->desc_, x.desc_));
  std::move(x.begin(), x.end(), begin());

  return *this;
}

template <typename T, std::size_t N>
MatrixRef<T, N> &MatrixRef<T, N>::operator=(const MatrixRef &x) {
  assert(same_extents(this->desc_, x.desc_));
  std::copy(x.begin(), x.end(), begin());

  return *this;
}

template <typename T, std::size_t N>
template <typename U>
MatrixRef<T, N>::MatrixRef(const MatrixRef<U, N> &x)
    : MatrixBase<T, N>{x.descriptor()}, ptr_(x.data()) {}

template <typename T, std::size_t N>
template <typename U>
MatrixRef<T, N> &MatrixRef<T, N>::operator=(const MatrixRef<U, N> &x) {
  static_assert(Convertible<U, T>(), "MatrixRef =: incompatible element types");
  assert(this->desc_.extents == x.descriptor().extents);

  std::copy(x.begin(), x.end(), begin());
  return *this;
}

template <typename T, std::size_t N>
template <typename U>
MatrixRef<T, N>::MatrixRef(const Matrix<U, N> &x)
    : MatrixBase<T, N>{x.descriptor()}, ptr_(x.data()) {}

template <typename T, std::size_t N>
template <typename U>
MatrixRef<T, N> &MatrixRef<T, N>::operator=(const Matrix<U, N> &x) {
  static_assert(Convertible<U, T>(), "MatrixRef =: incompatible element types");
  assert(this->desc_.extents == x.descriptor().extents);

  std::copy(x.begin(), x.end(), begin());
  return *this;
}

template <typename T, std::size_t N>
MatrixRef<T, N> &MatrixRef<T, N>::operator=(MatrixInitializer<T, N> init) {
  // std::array<std::size_t, N> extents = matrix_impl::derive_extents<N>(init);
  assert(matrix_impl::derive_extents<N>(init) == this->desc_.extents);

  auto iter = begin();
  matrix_impl::copy_flat(init, iter);

  return *this;
}

template <typename T, size_t N>
template <typename... Args>
Enable_if<matrix_impl::Requesting_slice<Args...>(), MatrixRef<T, N>>
MatrixRef<T, N>::operator()(const Args &... args) {
  MatrixSlice<N> d;
  d.start = this->desc_.start + matrix_impl::do_slice(this->desc_, d, args...);
  d.size = matrix_impl::compute_size(d.extents);
  return {d, data()};
}

template <typename T, size_t N>
template <typename... Args>
Enable_if<matrix_impl::Requesting_slice<Args...>(), MatrixRef<const T, N>>
MatrixRef<T, N>::operator()(const Args &... args) const {
  MatrixSlice<N> d;
  d.start = this->desc_.start + matrix_impl::do_slice(this->desc_, d, args...);
  d.size = matrix_impl::compute_size(d.extents);
  return {d, data()};
}

// row
template <typename T, size_t N>
MatrixRef<T, N - 1> MatrixRef<T, N>::row(size_t n) {
  assert(n < this->n_rows());
  MatrixSlice<N - 1> row;
  matrix_impl::slice_dim<0>(n, this->desc_, row);
  return {row, ptr_};
}

template <typename T, size_t N>
MatrixRef<const T, N - 1> MatrixRef<T, N>::row(size_t n) const {
  assert(n < this->n_rows());
  MatrixSlice<N - 1> row;
  matrix_impl::slice_dim<0>(n, this->desc_, row);
  return {row, ptr_};
}

// col
template <typename T, size_t N>
MatrixRef<T, N - 1> MatrixRef<T, N>::col(size_t n) {
  assert(n < this->cols());
  MatrixSlice<N - 1> col;
  matrix_impl::slice_dim<1>(n, this->desc_, col);
  return {col, ptr_};
}

template <typename T, size_t N>
MatrixRef<const T, N - 1> MatrixRef<T, N>::col(size_t n) const {
  assert(n < this->cols());
  MatrixSlice<N - 1> col;
  matrix_impl::slice_dim<1>(n, this->desc_, col);
  return {col, ptr_};
}

template <typename T, std::size_t N>
MatrixRef<T, N> MatrixRef<T, N>::rows(std::size_t i, std::size_t j) {
  assert(i < j);
  assert(j < this->n_rows());

  MatrixSlice<N> d;
  d.start = this->desc_.start;
  d.start += matrix_impl::do_slice_dim<N>(this->desc_, d, slice{i, j - i + 1});
  std::size_t NRest = N - 1;
  while (NRest >= 1) {
    d.start += matrix_impl::do_slice_dim2(this->desc_, d, slice{0}, NRest);
    --NRest;
  }
  return {d, data()};
}

template <typename T, std::size_t N>
MatrixRef<const T, N> MatrixRef<T, N>::rows(std::size_t i,
                                            std::size_t j) const {
  assert(i < j);
  assert(j < this->n_rows());

  MatrixSlice<N> d;
  d.start = this->desc_.start;
  d.start += matrix_impl::do_slice_dim<N>(this->desc_, d, slice{i, j - i + 1});
  std::size_t NRest = N - 1;
  while (NRest >= 1) {
    d.start += matrix_impl::do_slice_dim2(this->desc_, d, slice{0}, NRest);
    --NRest;
  }
  return {d, data()};
}

template <typename T, std::size_t N>
MatrixRef<T, N> MatrixRef<T, N>::cols(std::size_t i, std::size_t j) {
  assert(N >= 2);
  assert(i < j);
  assert(j < this->n_cols());

  MatrixSlice<N> d;
  d.start = this->desc_.start;
  d.start += matrix_impl::do_slice_dim<N>(this->desc_, d, slice{0});
  d.start +=
      matrix_impl::do_slice_dim<N - 1>(this->desc_, d, slice{i, j - i + 1});

  std::size_t NRest = N - 2;
  while (NRest >= 1) {
    d.start += matrix_impl::do_slice_dim2(this->desc_, d, slice{0}, NRest);
    --NRest;
  }
  return {d, data()};
}

template <typename T, std::size_t N>
MatrixRef<const T, N> MatrixRef<T, N>::cols(std::size_t i,
                                            std::size_t j) const {
  assert(N >= 2);
  assert(i < j);
  assert(j < this->n_cols());

  MatrixSlice<N> d;
  d.start = this->desc_.start;
  d.start += matrix_impl::do_slice_dim<N>(this->desc_, d, slice{0});
  d.start +=
      matrix_impl::do_slice_dim<N - 1>(this->desc_, d, slice{i, j - i + 1});

  std::size_t NRest = N - 2;
  while (NRest >= 1) {
    d.start += matrix_impl::do_slice_dim2(this->desc_, d, slice{0}, NRest);
    --NRest;
  }
  return {d, data()};
}

template <typename T, std::size_t N>
template <typename F>
MatrixRef<T, N> &MatrixRef<T, N>::apply(F f) {
  for (auto iter = begin(); iter != end(); ++iter) f(*iter);
  return *this;
}

template <typename T, std::size_t N>
template <typename M, typename F>
Enable_if<Matrix_type<M>(), MatrixRef<T, N> &> MatrixRef<T, N>::apply(
    const M &m, F f) {
  assert(same_extents(this->desc_, m.descriptor()));
  auto j = m.begin();
  for (auto i = begin(); i != end(); ++i) {
    f(*i, *j);
    ++j;
  }

  return *this;
}

template <typename T, std::size_t N>
MatrixRef<T, N> &MatrixRef<T, N>::operator=(const T &val) {
  return apply([&](T &a) { a = val; });
}

template <typename T, std::size_t N>
MatrixRef<T, N> &MatrixRef<T, N>::operator+=(const T &val) {
  return apply([&](T &a) { a += val; });
}

template <typename T, std::size_t N>
MatrixRef<T, N> &MatrixRef<T, N>::operator-=(const T &val) {
  return apply([&](T &a) { a -= val; });
}

template <typename T, std::size_t N>
MatrixRef<T, N> &MatrixRef<T, N>::operator*=(const T &val) {
  return apply([&](T &a) { a *= val; });
}

template <typename T, std::size_t N>
MatrixRef<T, N> &MatrixRef<T, N>::operator/=(const T &val) {
  return apply([&](T &a) { a /= val; });
}

template <typename T, std::size_t N>
MatrixRef<T, N> &MatrixRef<T, N>::operator%=(const T &val) {
  return apply([&](T &a) { a %= val; });
}

template <typename T, std::size_t N>
template <typename M>
Enable_if<Matrix_type<M>(), MatrixRef<T, N> &> MatrixRef<T, N>::operator+=(
    const M &m) {
  // static_assert(m.order_ == N, "+=: mismatched Matrix dimensions");
  assert(same_extents(this->desc_, m.descriptor()));  // make sure sizes match

  return apply(m, [&](T &a, const Value_type<M> &b) { a += b; });
}

template <typename T, std::size_t N>
template <typename M>
Enable_if<Matrix_type<M>(), MatrixRef<T, N> &> MatrixRef<T, N>::operator-=(
    const M &m) {
  // static_assert(m.order_ == N, "+=: mismatched Matrix dimensions");
  assert(same_extents(this->desc_, m.descriptor()));  // make sure sizes match

  return apply(m, [&](T &a, const Value_type<M> &b) { a -= b; });
}

template <typename T, std::size_t N>
template <typename M>
Enable_if<Matrix_type<M>(), MatrixRef<T, N> &> MatrixRef<T, N>::operator*=(
    const M &m) {
  assert(same_extents(this->desc_, m.descriptor()));  // make sure sizes match

  return apply(m, [&](T &a, const Value_type<M> &b) { a *= b; });
}

template <typename T, std::size_t N>
template <typename M>
Enable_if<Matrix_type<M>(), MatrixRef<T, N> &> MatrixRef<T, N>::operator/=(
    const M &m) {
  assert(same_extents(this->desc_, m.descriptor()));  // make sure sizes match

  return apply(m, [&](T &a, const Value_type<M> &b) { a /= b; });
}

template <typename T, std::size_t N>
template <typename M>
Enable_if<Matrix_type<M>(), MatrixRef<T, N> &> MatrixRef<T, N>::operator%=(
    const M &m) {
  assert(same_extents(this->desc_, m.descriptor()));  // make sure sizes match

  return apply(m, [&](T &a, const Value_type<M> &b) { a %= b; });
}

template <typename T, std::size_t N>
template <typename U>
Matrix<U, N> MatrixRef<T, N>::operator-() const {
  Matrix<U, N> res(*this);
  return res.apply([&](U &a) { a = -a; });
}

 template <typename T, std::size_t N>
 template <typename U, std::size_t NN, typename X>
   MatrixRef<T, N> &MatrixRef<T, N>::operator=(const MatrixRef<U, 2> &x) {
   static_assert(Convertible<U, T>(),
                 "MatrixRef =: incompatible element types");
   assert(this->size() == x.size());
   assert(x.n_cols() == 1);

   std::copy(x.begin(), x.end(), begin());

   return *this;
 }

 template <typename T, std::size_t N>
 template <typename U, std::size_t NN, typename X>
   MatrixRef<T, N> &MatrixRef<T, N>::operator=(const Matrix<U, 2> &x) {
   static_assert(Convertible<U, T>(),
                 "MatrixRef =: incompatible element types");
   assert(this->size() == x.size());
   assert(x.n_cols() == 1);

   std::copy(x.begin(), x.end(), begin());

   return *this;
 }

 template <typename T, std::size_t N>
 template <typename U, std::size_t NN, typename X>
   MatrixRef<T, N> &MatrixRef<T, N>::operator=(const MatrixRef<U, 1> &x) {
   static_assert(Convertible<U, T>(),
                 "MatrixRef =: incompatible element types");
   assert(this->size() == x.size());
   assert(this->n_cols() == 1);

   std::copy(x.begin(), x.end(), begin());

   return *this;
 }

 template <typename T, std::size_t N>
   template <typename U, std::size_t NN, typename X>
   MatrixRef<T, N> &MatrixRef<T, N>::operator=(const Matrix<U, 1> &x) {
   static_assert(Convertible<U, T>(),
                 "MatrixRef =: incompatible element types");
   assert(this->size() == x.size());
   assert(this->n_cols() == 1);

   std::copy(x.begin(), x.end(), begin());

   return *this;
 }

template <typename T>
class MatrixRef<T, 0> : public MatrixBase<T, 0> {
 public:
  using iterator = T *;
  using const_iterator = const T *;

  MatrixRef(const MatrixSlice<0> &s, T *p) : ptr_{p + s.start} {}

  //! total number of elements
  std::size_t size() const { return 1; }

  //! "flat" element access
  ///@{
  T *data() { return ptr_; }
  const T *data() const { return ptr_; }
  ///@}

  T &operator()() { return *ptr_; };
  const T &operator()() const { return *ptr_; }

  operator T &() { return *ptr_; }
  operator const T &() const { return *ptr_; }

  iterator begin() { return iterator(data()); }
  const_iterator begin() const { return const_iterator(data()); }
  iterator end() { return iterator(data() + 1); }
  const_iterator end() const { return const_iterator(data() + 1); }

 private:
  T *ptr_;
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const MatrixRef<T, 0> &mr0) {
  return os << (const T &)mr0;
}

template <typename T, std::size_t N>
class MatrixRefIterator {
  template <typename U, size_t NN>
  friend std::ostream &operator<<(std::ostream &os,
                                  const MatrixRefIterator<U, NN> &iter);

 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = typename std::remove_const<T>::type;
  using pointer = T *;
  using reference = T &;
  using difference_type = std::ptrdiff_t;

  MatrixRefIterator(const MatrixSlice<N> &s, T *base, bool limit = false);
  MatrixRefIterator &operator=(const MatrixRefIterator &);

  const MatrixSlice<N> &descriptor() const { return desc_; }

  T &operator*() { return *ptr_; }
  T *operator->() { return ptr_; }

  const T &operator*() const { return *ptr_; }
  const T *operator->() const { return ptr_; }

  MatrixRefIterator &operator++();
  MatrixRefIterator operator++(int);

 private:
  void increment();

  std::array<size_t, N> indx_;
  const MatrixSlice<N> &desc_;
  T *ptr_;
};

template <typename T, std::size_t N>
MatrixRefIterator<T, N>::MatrixRefIterator(const MatrixSlice<N> &s, T *base,
                                           bool limit)
    : desc_(s) {
  std::fill(indx_.begin(), indx_.end(), 0);

  if (limit) {
    indx_[0] = desc_.extents[0];
    ptr_ = base + desc_.offset(indx_);
  } else {
    ptr_ = base + s.start;
  }
}

template <typename T, std::size_t N>
MatrixRefIterator<T, N> &MatrixRefIterator<T, N>::operator=(
    const MatrixRefIterator &iter) {
  std::copy(iter.indx_.begin(), iter.indx_.end(), indx_.begin());
  ptr_ = iter.ptr_;

  return *this;
}

template <typename T, std::size_t N>
MatrixRefIterator<T, N> &MatrixRefIterator<T, N>::operator++() {
  increment();
  return *this;
}

template <typename T, std::size_t N>
MatrixRefIterator<T, N> MatrixRefIterator<T, N>::operator++(int) {
  MatrixRefIterator<T, N> x = *this;
  increment();
  return *x;
}

template <typename T, std::size_t N>
void MatrixRefIterator<T, N>::increment() {
  std::size_t d = N - 1;

  while (true) {
    ptr_ += desc_.strides[d];
    ++indx_[d];

    if (indx_[d] != desc_.extents[d]) break;

    if (d != 0) {
      ptr_ -= desc_.strides[d] * desc_.extents[d];
      indx_[d] = 0;
      --d;
    } else {
      break;
    }
  }
}

template <typename T, size_t N>
std::ostream &operator<<(std::ostream &os,
                         const MatrixRefIterator<T, N> &iter) {
  os << "target: " << *iter.ptr_ << ", indx: " << iter.indx_ << std::endl;
  return os;
}

template <typename T, std::size_t N>
inline bool operator==(const MatrixRefIterator<T, N> &a,
                       const MatrixRefIterator<T, N> &b) {
  assert(a.descriptor() == b.descriptor());
  return &*a == &*b;
}

template <typename T, std::size_t N>
inline bool operator!=(const MatrixRefIterator<T, N> &a,
                       const MatrixRefIterator<T, N> &b) {
  return !(a == b);
}

}  // namespace slab

#endif  // SLAB_MATRIX_MATRIX_REF_H_
