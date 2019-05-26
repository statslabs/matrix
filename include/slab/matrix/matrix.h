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

/// @file matrix.h
/// @brief A Matrix template

#ifndef _SLAB_MATRIX_MATRIX_H
#define _SLAB_MATRIX_MATRIX_H

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include "slab/__config"

#include "slab/matrix/matrix_base.h"
#include "slab/matrix/matrix_fwd.h"
#include "slab/matrix/matrix_ref.h"
#include "slab/matrix/matrix_slice.h"
#include "slab/matrix/packed_matrix.h"
#include "slab/matrix/support.h"

_SLAB_BEGIN_NAMESPACE

//! Matrix<T,N> is an N-dimensional matrix of some value type T.
/*!
 * \tparam T value type.
 * \tparam N number of dimensions.
 *
 * This class implements matrix class which provides support subscripting,
 * slicing and basic matrix arithmetic operations.
 */
template <typename T, std::size_t N>
class Matrix : public MatrixBase<T, N> {
  // ----------------------------------------
  // The core member functions in book 'TCPL'
  // ----------------------------------------
 public:
  //! @cond Doxygen_Suppress
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;

  Matrix() = default;
  Matrix(Matrix &&) = default;  // move
  Matrix &operator=(Matrix &&) = default;
  Matrix(const Matrix &) = default;  // copy
  Matrix &operator=(const Matrix &) = default;
  ~Matrix() = default;
  //! @endcond

  //! construct from Matrix
  template <typename M, typename = Enable_if<Matrix_type<M>()>>
  Matrix(const M &x);
  //! assign from Matrix
  template <typename M, typename = Enable_if<Matrix_type<M>()>>
  Matrix &operator=(const M &x);

  //! construct from MatrixRef
  template <typename U>
  Matrix(const MatrixRef<U, N> &);
  //! assign from MatrixRef
  template <typename U>
  Matrix &operator=(const MatrixRef<U, N> &);

  //! specify the extents
  template <typename... Exts>
  explicit Matrix(Exts... exts);

  //! initialize from list
  Matrix(MatrixInitializer<T, N>);
  //! assign from list
  Matrix &operator=(MatrixInitializer<T, N>);

  //! don't use {} except for elements
  ///@{
  template <typename U, std::size_t NN = N,
            typename = Enable_if<Convertible<U, std::size_t>()>,
            typename = Enable_if<(NN > 1)>>
  Matrix(std::initializer_list<U>) = delete;
  template <typename U, std::size_t NN = N,
            typename = Enable_if<Convertible<U, std::size_t>()>,
            typename = Enable_if<(NN > 1)>>
  Matrix &operator=(std::initializer_list<U>) = delete;
  ///@}

  //! total number of elements
  std::size_t size() const { return elems_.size(); }

  //! "flat" element access
  ///@{
  T *data() { return elems_.data(); }
  const T *data() const { return elems_.data(); }
  ///@}

 private:
  std::vector<T> elems_;  // the elements

  // ---------------------------------------------
  // Member functions for subscripting and slicing
  // ---------------------------------------------
 public:
  //! m(i,j,k) subscripting with integers
  using MatrixBase<T, N>::operator();

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
  iterator begin() { return elems_.begin(); }
  const_iterator begin() const { return elems_.cbegin(); }
  iterator end() { return elems_.end(); }
  const_iterator end() const { return elems_.cend(); }
  ///@}

  // --------------------------------------------------
  // Member functions for matrix arithmetic operations
  // --------------------------------------------------
 public:
  //! matrix arithmetic operations
  ///@{
  template <typename F>
  Matrix &apply(F f);  // f(x) for every element x

  // f(x, mx) for corresponding elements of *this and m
  template <typename M, typename F>
  Enable_if<Matrix_type<M>(), Matrix &> apply(const M &m, F f);

  Matrix &operator=(const T &value);  // assignment with scalar

  Matrix &operator+=(const T &value);  // scalar addition
  Matrix &operator-=(const T &value);  // scalar subtraction
  Matrix &operator*=(const T &value);  // scalar multiplication
  Matrix &operator/=(const T &value);  // scalar division
  Matrix &operator%=(const T &value);  // scalar modulo

  // matrix addition
  template <typename M>
  Enable_if<Matrix_type<M>(), Matrix &> operator+=(const M &x);
  // matrix subtraction
  template <typename M>
  Enable_if<Matrix_type<M>(), Matrix &> operator-=(const M &x);
  // element-wise multiplication
  template <typename M>
  Enable_if<Matrix_type<M>(), Matrix &> operator*=(const M &x);
  // element-wise division
  template <typename M>
  Enable_if<Matrix_type<M>(), Matrix &> operator/=(const M &x);
  // element-wise modulus
  template <typename M>
  Enable_if<Matrix_type<M>(), Matrix &> operator%=(const M &x);

  template <typename U = typename std::remove_const<T>::type>
  Matrix<U, N> operator-() const;
  ///@}

  // -----------------------------------
  // Member functions for a Matrix<T, 1>
  // -----------------------------------

 public:
  //! Construct a vector from a matrix
  ///@{
  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 1)>>
  Matrix(const Matrix<U, 2> &x);
  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 1)>>
  Matrix(const MatrixRef<U, 2> &x);
  ///@}

  //！Assign a vector from a matrix
  ///@{
  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 1)>>
  Matrix &operator=(const Matrix<U, 2> &x);
  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 1)>>
  Matrix &operator=(const MatrixRef<U, 2> &x);
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

  //! Construct a std::vector from a slab::vec
  template <std::size_t NN = N, typename = Enable_if<(NN == 1)>>
  std::vector<T> std_vec() const {
    return std::vector<T>(begin(), end());
  }

  // -----------------------------------
  // Member functions for a Matrix<T, 2>
  // -----------------------------------

 public:
  //! Construct a matrix from a vector
  ///@{
  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 2)>>
  Matrix(const Matrix<U, 1> &x);
  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 2)>>
  Matrix(const MatrixRef<U, 1> &x);
  ///@}

  //! Assign a matrix from a vector
  ///@{
  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 2)>>
  Matrix &operator=(const Matrix<U, 1> &x);
  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 2)>>
  Matrix &operator=(const MatrixRef<U, 1> &x);
  ///@}

  //! Construct a matrix from a symmetric/triangular/hermitian matrix
  ///@{
  template <typename U, typename TRI, std::size_t NN = N,
            typename = Enable_if<(NN == 2)>>
  Matrix(const SymmetricMatrix<U, TRI> &x);

  template <typename U, typename TRI, std::size_t NN = N,
            typename = Enable_if<(NN == 2)>>
  Matrix(const TriangularMatrix<U, TRI> &x);

  template <typename U, typename TRI, std::size_t NN = N,
            typename = Enable_if<(NN == 2)>>
  Matrix(const HermitianMatrix<U, TRI> &x);
  ///@}

  //! Assign a matrix from a symmetric/triangular/hermitian matrix
  ///@{
  template <typename U, typename TRI, std::size_t NN = N,
            typename = Enable_if<(NN == 2)>>
  Matrix &operator=(const SymmetricMatrix<U, TRI> &x);

  template <typename U, typename TRI, std::size_t NN = N,
            typename = Enable_if<(NN == 2)>>
  Matrix &operator=(const TriangularMatrix<U, TRI> &x);

  template <typename U, typename TRI, std::size_t NN = N,
            typename = Enable_if<(NN == 2)>>
  Matrix &operator=(const HermitianMatrix<U, TRI> &x);
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

  //! print mat/vec object to std::cout
  template <std::size_t NN = N, typename = Enable_if<(NN == 1) || (NN == 2)>>
  void print(const std::string &str = "") const {
    printf("\n %s\n", str.c_str());
    for (std::size_t i = 0; i != this->n_rows(); ++i) {
      raw_print(row(i));
      printf("\n");
    }
  }

  // ---------------------------------------------
  // More member functions for the matrix template
  // ---------------------------------------------

 public:
  Matrix(const Matrix<double, N> &x, const Matrix<double, N> &y)
      : MatrixBase<T, N>(x.descriptor()) {
    _SLAB_ASSERT(is_complex_double<T>::value, "invalid constructor");
    _SLAB_ASSERT(x.descriptor() == y.descriptor(), "x and y size unmatched");

    for (auto iter1 = x.begin(), iter2 = y.begin(); iter1 != x.end();
         ++iter1, ++iter2) {
      elems_.push_back(std::complex<double>(*iter1, *iter2));
    }
  }

  Matrix(const Matrix<float, N> &x, const Matrix<float, N> &y)
      : MatrixBase<T, N>(x.descriptor()) {
    _SLAB_ASSERT(is_complex_float<T>::value, "invalid constructor");
    _SLAB_ASSERT(x.descriptor() == y.descriptor(), "x and y size unmatched");

    for (auto iter1 = x.begin(), iter2 = y.begin(); iter1 != x.end();
         ++iter1, ++iter2) {
      elems_.push_back(std::complex<float>(*iter1, *iter2));
    }
  }

  //! clear content
  void clear();

  //! return matrix transpose
  template <std::size_t NN = N, typename = Enable_if<(NN == 1) || (NN == 2)>>
  Matrix<T, 2> t() const {
    return transpose(*this);
  }

  //! return inverse of square matrix
  template <std::size_t NN = N, typename = Enable_if<(NN == 2)>>
  Matrix<T, 2> i() const {
    return inverse(*this);
  }

  //! check whether object is empty
  ///@{
  bool empty() const { return begin() == end(); }
  bool is_empty() const { return empty(); }
  ///@}

  //  void save(const std::string &filename) {
  //    std::ostream os(filename);
  //  }
  //! load matrix from a file
  void load(const std::string &filename);

#ifdef _SLAB_USE_RCPP_AS_WRAP
  // -----------------------------
  // Conversion between R and C++
  // -----------------------------
 public:
  // this ctor enables implicit Rcpp::as
  Matrix(SEXP s);

  // this operator enables implicit Rcpp::wrap
  operator SEXP();
  operator SEXP() const;

#endif
};

//! @cond Doxygen_Suppress

template <typename T, std::size_t N>
template <typename M, typename X>
Matrix<T, N>::Matrix(const M &x)
    : MatrixBase<T, N>(x.descriptor()), elems_(x.begin(), x.end()) {
  static_assert(Convertible<typename M::value_type, T>(), "");
}

template <typename T, std::size_t N>
template <typename M, typename X>
Matrix<T, N> &Matrix<T, N>::operator=(const M &x) {
  static_assert(Convertible<typename M::value_type, T>(), "");

  this->desc_ = x.descriptor();
  elems_.assign(x.begin(), x.end());
  return *this;
}

template <typename T, std::size_t N>
template <typename U>
Matrix<T, N>::Matrix(const MatrixRef<U, N> &x)  // copy desc_ and elements
    : MatrixBase<T, N>{x.descriptor().extents}, elems_{x.begin(), x.end()} {
  static_assert(Convertible<U, T>(),
                "Matrix constructor: incompatible element types");
}

template <typename T, std::size_t N>
template <typename U>
Matrix<T, N> &Matrix<T, N>::operator=(const MatrixRef<U, N> &x) {
  static_assert(Convertible<U, T>(), "Matrix =: incompatible element types");

  this->desc_ = x.descriptor();
  this->desc_.start = 0;
  elems_.assign(x.begin(), x.end());
  return *this;
}

template <typename T, std::size_t N>
template <typename... Exts>
Matrix<T, N>::Matrix(Exts... exts)
    : MatrixBase<T, N>{exts...},  // copy extents
      elems_(this->desc_.size)    // allocate desc_.size elements and initialize
{}

template <typename T, std::size_t N>
Matrix<T, N>::Matrix(MatrixInitializer<T, N> init) {
  // intialize start
  this->desc_.start = 0;
  // deduce extents from initializer list
  this->desc_.extents = matrix_impl::derive_extents<N>(init);
  // compute strides and size
  this->desc_.size =
      matrix_impl::compute_strides(this->desc_.extents, this->desc_.strides);

  elems_.reserve(this->desc_.size);        // make room for slices
  matrix_impl::insert_flat(init, elems_);  // initialize from initializer list
  assert(elems_.size() == this->desc_.size);
}

template <typename T, std::size_t N>
Matrix<T, N> &Matrix<T, N>::operator=(MatrixInitializer<T, N> init) {
  elems_.clear();

  // intialize start
  this->desc_.start = 0;
  // deduce extents from initializer list
  this->desc_.extents = matrix_impl::derive_extents<N>(init);
  // compute strides and size
  this->desc_.size =
      matrix_impl::compute_strides(this->desc_.extents, this->desc_.strides);

  elems_.reserve(this->desc_.size);        // make room for slices
  matrix_impl::insert_flat(init, elems_);  // initialize from initializer list
  assert(elems_.size() == this->desc_.size);

  return *this;
}

template <typename T, std::size_t N>
template <typename... Args>
Enable_if<matrix_impl::Requesting_slice<Args...>(), MatrixRef<T, N>>
Matrix<T, N>::operator()(const Args &... args) {
  MatrixSlice<N> d;
  d.start = matrix_impl::do_slice(this->desc_, d, args...);
  d.size = matrix_impl::compute_size(d.extents);
  return {d, data()};
}

template <typename T, std::size_t N>
template <typename... Args>
Enable_if<matrix_impl::Requesting_slice<Args...>(), MatrixRef<const T, N>>
Matrix<T, N>::operator()(const Args &... args) const {
  MatrixSlice<N> d;
  d.start = matrix_impl::do_slice(this->desc_, d, args...);
  d.size = matrix_impl::compute_size(d.extents);
  return {d, data()};
}

// row
template <typename T, std::size_t N>
MatrixRef<T, N - 1> Matrix<T, N>::row(std::size_t n) {
  assert(n < this->n_rows());
  MatrixSlice<N - 1> row;
  matrix_impl::slice_dim<0>(n, this->desc_, row);
  return {row, data()};
}

template <typename T, std::size_t N>
MatrixRef<const T, N - 1> Matrix<T, N>::row(std::size_t n) const {
  assert(n < this->n_rows());
  MatrixSlice<N - 1> row;
  matrix_impl::slice_dim<0>(n, this->desc_, row);
  return {row, data()};
}

// col
template <typename T, std::size_t N>
MatrixRef<T, N - 1> Matrix<T, N>::col(std::size_t n) {
  assert(n < this->n_cols());
  MatrixSlice<N - 1> col;
  matrix_impl::slice_dim<1>(n, this->desc_, col);
  return {col, data()};
}

template <typename T, std::size_t N>
MatrixRef<const T, N - 1> Matrix<T, N>::col(std::size_t n) const {
  assert(n < this->n_cols());
  MatrixSlice<N - 1> col;
  matrix_impl::slice_dim<1>(n, this->desc_, col);
  return {col, data()};
}

template <typename T, std::size_t N>
MatrixRef<T, N> Matrix<T, N>::rows(std::size_t i, std::size_t j) {
  assert(i <= j);
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
MatrixRef<const T, N> Matrix<T, N>::rows(std::size_t i, std::size_t j) const {
  assert(i <= j);
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
MatrixRef<T, N> Matrix<T, N>::cols(std::size_t i, std::size_t j) {
  assert(N >= 2);
  assert(i <= j);
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
MatrixRef<const T, N> Matrix<T, N>::cols(std::size_t i, std::size_t j) const {
  assert(N >= 2);
  assert(i <= j);
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
Matrix<T, N> &Matrix<T, N>::apply(F f) {
  for (auto &x : elems_) f(x);  // this loop uses stride iterators
  return *this;
}

template <typename T, std::size_t N>
template <typename M, typename F>
Enable_if<Matrix_type<M>(), Matrix<T, N> &> Matrix<T, N>::apply(const M &m,
                                                                F f) {
  assert(same_extents(this->desc_, m.descriptor()));
  auto j = m.begin();
  for (auto i = begin(); i != end(); ++i) {
    f(*i, *j);
    ++j;
  }

  return *this;
}

template <typename T, std::size_t N>
Matrix<T, N> &Matrix<T, N>::operator=(const T &val) {
  return apply([&](T &a) { a = val; });
}

template <typename T, std::size_t N>
Matrix<T, N> &Matrix<T, N>::operator+=(const T &val) {
  return apply([&](T &a) { a += val; });
}

template <typename T, std::size_t N>
Matrix<T, N> &Matrix<T, N>::operator-=(const T &val) {
  return apply([&](T &a) { a -= val; });
}

template <typename T, std::size_t N>
Matrix<T, N> &Matrix<T, N>::operator*=(const T &val) {
  return apply([&](T &a) { a *= val; });
}

template <typename T, std::size_t N>
Matrix<T, N> &Matrix<T, N>::operator/=(const T &val) {
  return apply([&](T &a) { a /= val; });
}

template <typename T, std::size_t N>
Matrix<T, N> &Matrix<T, N>::operator%=(const T &val) {
  return apply([&](T &a) { a %= val; });
}

template <typename T, std::size_t N>
template <typename M>
Enable_if<Matrix_type<M>(), Matrix<T, N> &> Matrix<T, N>::operator+=(
    const M &m) {
  // static_assert(m.order_ == N, "+=: mismatched Matrix dimensions");
  assert(same_extents(this->desc_, m.descriptor()));  // make sure sizes match

  return apply(m, [&](T &a, const Value_type<M> &b) { a += b; });
}

template <typename T, std::size_t N>
template <typename M>
Enable_if<Matrix_type<M>(), Matrix<T, N> &> Matrix<T, N>::operator-=(
    const M &m) {
  // static_assert(m.order_ == N, "-=: mismatched Matrix dimensions");
  assert(same_extents(this->desc_, m.descriptor()));  // make sure sizes match

  return apply(m, [&](T &a, const Value_type<M> &b) { a -= b; });
}

template <typename T, std::size_t N>
template <typename M>
Enable_if<Matrix_type<M>(), Matrix<T, N> &> Matrix<T, N>::operator*=(
    const M &m) {
  assert(same_extents(this->desc_, m.descriptor()));  // make sure sizes match

  return apply(m, [&](T &a, const Value_type<M> &b) { a *= b; });
}

template <typename T, std::size_t N>
template <typename M>
Enable_if<Matrix_type<M>(), Matrix<T, N> &> Matrix<T, N>::operator/=(
    const M &m) {
  assert(same_extents(this->desc_, m.descriptor()));  // make sure sizes match

  return apply(m, [&](T &a, const Value_type<M> &b) { a /= b; });
}

template <typename T, std::size_t N>
template <typename M>
Enable_if<Matrix_type<M>(), Matrix<T, N> &> Matrix<T, N>::operator%=(
    const M &m) {
  assert(same_extents(this->desc_, m.descriptor()));  // make sure sizes match

  return apply(m, [&](T &a, const Value_type<M> &b) { a %= b; });
}

template <typename T, std::size_t N>
template <typename U>
Matrix<U, N> Matrix<T, N>::operator-() const {
  Matrix<U, N> res(*this);
  return res.apply([&](T &a) { a = -a; });
}

template <typename T, std::size_t N>
template <typename U, std::size_t NN, typename X>
Matrix<T, N>::Matrix(const Matrix<U, 2> &x)
    : MatrixBase<T, N>{x.n_rows()}, elems_{x.begin(), x.end()} {
  static_assert(Convertible<U, T>(),
                "Matrix constructor: incompatible element types");
  assert(x.n_cols() == 1);
}

template <typename T, std::size_t N>
template <typename U, std::size_t NN, typename X>
Matrix<T, N>::Matrix(const MatrixRef<U, 2> &x)
    : MatrixBase<T, N>{x.n_rows()}, elems_{x.begin(), x.end()} {
  static_assert(Convertible<U, T>(),
                "Matrix constructor: incompatible element types");
  assert(x.n_cols() == 1);
}

template <typename T, std::size_t N>
template <typename U, std::size_t NN, typename X>
Matrix<T, N> &Matrix<T, N>::operator=(const Matrix<U, 2> &x) {
  static_assert(Convertible<U, T>(), "Matrix =: incompatible element types");
  assert(x.n_cols() == 1);

  this->desc_.size = x.descriptor().size;
  this->desc_.start = 0;
  this->desc_.extents[0] = x.n_rows();
  this->desc_.strides[0] = 1;

  elems_.assign(x.begin(), x.end());

  return *this;
}

template <typename T, std::size_t N>
template <typename U, std::size_t NN, typename X>
Matrix<T, N> &Matrix<T, N>::operator=(const MatrixRef<U, 2> &x) {
  static_assert(Convertible<U, T>(), "Matrix =: incompatible element types");
  assert(x.n_cols() == 1);

  this->desc_.size = x.descriptor().size;
  this->desc_.start = 0;
  this->desc_.extents[0] = x.n_rows();
  this->desc_.strides[0] = 1;

  elems_.assign(x.begin(), x.end());

  return *this;
}

template <typename T, std::size_t N>
template <typename U, std::size_t NN, typename X>
Matrix<T, N>::Matrix(const Matrix<U, 1> &x)
    : MatrixBase<T, N>{x.n_rows(), 1}, elems_{x.begin(), x.end()} {
  static_assert(Convertible<U, T>(),
                "Matrix constructor: incompatible element types");
}

template <typename T, std::size_t N>
template <typename U, std::size_t NN, typename X>
Matrix<T, N>::Matrix(const MatrixRef<U, 1> &x)
    : MatrixBase<T, N>{x.n_rows(), 1}, elems_{x.begin(), x.end()} {
  static_assert(Convertible<U, T>(),
                "Matrix constructor: incompatible element types");
}

template <typename T, std::size_t N>
template <typename U, std::size_t NN, typename X>
Matrix<T, N> &Matrix<T, N>::operator=(const Matrix<U, 1> &x) {
  static_assert(Convertible<U, T>(), "Matrix =: incompatible element types");

  this->desc_.size = x.descriptor().size;
  this->desc_.start = 0;
  this->desc_.extents[0] = x.n_rows();
  this->desc_.extents[1] = 1;
  this->desc_.strides[0] = x.n_rows();
  this->desc_.strides[1] = 1;

  elems_.assign(x.begin(), x.end());

  return *this;
}

template <typename T, std::size_t N>
template <typename U, std::size_t NN, typename X>
Matrix<T, N> &Matrix<T, N>::operator=(const MatrixRef<U, 1> &x) {
  static_assert(Convertible<U, T>(), "Matrix =: incompatible element types");

  this->desc_.size = x.descriptor().size;
  this->desc_.start = 0;
  this->desc_.extents[0] = x.n_rows();
  this->desc_.extents[1] = 1;
  this->desc_.strides[0] = x.n_rows();
  this->desc_.strides[1] = 1;

  elems_.assign(x.begin(), x.end());

  return *this;
}

template <typename T, std::size_t N>
template <typename U, typename TRI, std::size_t NN, typename X>
Matrix<T, N>::Matrix(const SymmetricMatrix<U, TRI> &x)
    : MatrixBase<T, N>{x.n_rows(), x.n_cols()} {
  static_assert(Convertible<U, T>(),
                "Matrix constructor: incompatible element types");
  for (std::size_t i = 0; i != x.n_rows(); ++i) {
    for (std::size_t j = 0; j != x.n_cols(); ++j) {
      elems_.push_back(x(i, j));
    }
  }
}

template <typename T, std::size_t N>
template <typename U, typename TRI, std::size_t NN, typename X>
Matrix<T, N>::Matrix(const TriangularMatrix<U, TRI> &x)
    : MatrixBase<T, N>{x.n_rows(), x.n_cols()} {
  static_assert(Convertible<U, T>(),
                "Matrix constructor: incompatible element types");
  for (std::size_t i = 0; i != x.n_rows(); ++i) {
    for (std::size_t j = 0; j != x.n_cols(); ++j) {
      elems_.push_back(x(i, j));
    }
  }
}

template <typename T, std::size_t N>
template <typename U, typename TRI, std::size_t NN, typename X>
Matrix<T, N>::Matrix(const HermitianMatrix<U, TRI> &x)
    : MatrixBase<T, N>{x.n_rows(), x.n_cols()} {
  static_assert(Convertible<U, T>(),
                "Matrix constructor: incompatible element types");
  for (std::size_t i = 0; i != x.n_rows(); ++i) {
    for (std::size_t j = 0; j != x.n_cols(); ++j) {
      elems_.push_back(x(i, j));
    }
  }
}

template <typename T, std::size_t N>
template <typename U, typename TRI, std::size_t NN, typename X>
Matrix<T, N> &Matrix<T, N>::operator=(const SymmetricMatrix<U, TRI> &x) {
  static_assert(Convertible<U, T>(), "Matrix =: incompatible element types");

  std::size_t n = x.n_rows();

  this->desc_.size = n * n;
  this->desc_.start = 0;
  this->desc_.extents[0] = n;
  this->desc_.extents[1] = n;
  this->desc_.strides[0] = n;
  this->desc_.strides[1] = 1;

  elems_.reserve(n * n);

  for (std::size_t i = 0; i != x.n_rows(); ++i) {
    for (std::size_t j = 0; j != x.n_cols(); ++j) {
      *(data() + this->desc_(i, j)) = x(i, j);
    }
  }

  return *this;
}

template <typename T, std::size_t N>
template <typename U, typename TRI, std::size_t NN, typename X>
Matrix<T, N> &Matrix<T, N>::operator=(const TriangularMatrix<U, TRI> &x) {
  static_assert(Convertible<U, T>(), "Matrix =: incompatible element types");

  std::size_t n = x.n_rows();

  this->desc_.size = n * n;
  this->desc_.start = 0;
  this->desc_.extents[0] = n;
  this->desc_.extents[1] = n;
  this->desc_.strides[0] = n;
  this->desc_.strides[1] = 1;

  elems_.reserve(n * n);

  for (std::size_t i = 0; i != x.n_rows(); ++i) {
    for (std::size_t j = 0; j != x.n_cols(); ++j) {
      *(data() + this->desc_(i, j)) = x(i, j);
    }
  }

  return *this;
}

template <typename T, std::size_t N>
template <typename U, typename TRI, std::size_t NN, typename X>
Matrix<T, N> &Matrix<T, N>::operator=(const HermitianMatrix<U, TRI> &x) {
  static_assert(Convertible<U, T>(), "Matrix =: incompatible element types");

  std::size_t n = x.n_rows();

  this->desc_.size = n * n;
  this->desc_.start = 0;
  this->desc_.extents[0] = n;
  this->desc_.extents[1] = n;
  this->desc_.strides[0] = n;
  this->desc_.strides[1] = 1;

  elems_.reserve(n * n);

  for (std::size_t i = 0; i != x.n_rows(); ++i) {
    for (std::size_t j = 0; j != x.n_cols(); ++j) {
      *(data() + this->desc_(i, j)) = x(i, j);
    }
  }

  return *this;
}

template <typename T, std::size_t N>
void Matrix<T, N>::clear() {
  this->desc_.clear();
  elems_.clear();
}

template <typename T, std::size_t N>
void Matrix<T, N>::load(const std::string &filename) {
  std::ifstream is(filename);
  if (is.is_open()) {
    // read the first line
    std::string first_line;
    getline(is, first_line);

    // read the extents into ivec
    std::istringstream iss(first_line);
    int val;
    std::vector<int> ivec;
    while (iss >> val) ivec.push_back(val);

    if (ivec.size() != this->order())
      std::cout << "incorrect extents" << std::endl;
    this->desc_.start = 0;
    std::copy(ivec.begin(), ivec.end(), this->desc_.extents.begin());
    this->desc_.size =
        matrix_impl::compute_strides(this->desc_.extents, this->desc_.strides);

    std::istream_iterator<T> in(is), end;
    elems_.assign(in, end);
  } else {
    std::cout << "Fail to open the file" << std::endl;
  }
}

//! @endcond

#ifdef _SLAB_USE_RCPP_AS_WRAP

template <typename T, std::size_t N>
Matrix<T, N>::Matrix(SEXP s) {
  SEXP dims = Rf_getAttrib(s, R_DimSymbol);
  _SLAB_ASSERT(Rf_length(dims) == this->order(),
               "Matrix(SEXP): unmatched dimensions.");

  int num_dims = Rf_length(dims);
  int num_elems = Rf_length(s);
  SEXP s2 = PROTECT(Rf_allocVector(TYPEOF(s), (R_xlen_t)num_elems));
  Rf_setAttrib(s2, R_DimSymbol, dims);
  if (num_dims <= 2) Rf_copyMatrix(s2, s, TRUE);

  elems_.reserve(num_elems);
  for (int i = 0; i != num_elems; ++i) {
    if (Rf_isReal(s))
      elems_.push_back(REAL(s2)[i]);
    else if (Rf_isInteger(s))
      elems_.push_back(INTEGER(s2)[i]);
    else
      _SLAB_ERROR("Matrix(SEXP): unsupported SEXP type");
  }

  std::array<std::size_t, N> exts;
  if (Rf_isNull(dims))
    exts[0] = num_elems;
  else {
    exts[N - 2] = INTEGER(dims)[0];
    exts[N - 1] = INTEGER(dims)[1];
    for (int i = 2; i != N; ++i) exts[i] = INTEGER(dims)[i];
  }
  this->desc_ = MatrixSlice<N>(exts);

  UNPROTECT(1);
}

template <typename T, std::size_t N>
Matrix<T, N>::operator SEXP() {
  int num_elems = this->size();
  SEXP res = PROTECT(Rcpp::wrap(this->data(), this->data() + num_elems));
  SEXP res2 = PROTECT(Rf_allocVector(Rcpp::traits::r_sexptype_traits<T>::rtype,
                                     (R_xlen_t)num_elems));

  int num_dims = this->order();
  SEXP dim = PROTECT(Rf_allocVector(INTSXP, num_dims));
  int *idim = INTEGER(dim);
  if (num_dims <= 2) {
    for (int i = 0; i != num_dims; ++i) idim[i] = this->desc_.extents[i];
  } else {
    idim[0] = this->desc_.extents[num_dims - 2];
    idim[1] = this->desc_.extents[num_dims - 1];
    for (int i = 2; i != num_dims; ++i)
      idim[i] = this->desc_.extents[num_dims - i - 1];
  }
  Rf_setAttrib(res, R_DimSymbol, dim);
  Rf_setAttrib(res2, R_DimSymbol, dim);

  if (num_dims <= 2) Rf_copyMatrix(res2, res, TRUE);
  // TODO: else CopyMatrixByRow(res2, res); // for num_dims > 2

  UNPROTECT(3);

  return res2;
}

template <typename T, std::size_t N>
Matrix<T, N>::operator SEXP() const {
  int num_elems = this->size();
  SEXP res = PROTECT(Rcpp::wrap(this->data(), this->data() + num_elems));
  SEXP res2 = PROTECT(Rf_allocVector(Rcpp::traits::r_sexptype_traits<T>::rtype,
                                     (R_xlen_t)num_elems));

  int num_dims = this->order();
  SEXP dim = PROTECT(Rf_allocVector(INTSXP, num_dims));
  int *idim = INTEGER(dim);
  if (num_dims <= 2) {
    for (int i = 0; i != num_dims; ++i) idim[i] = this->desc_.extents[i];
  } else {
    idim[0] = this->desc_.extents[num_dims - 2];
    idim[1] = this->desc_.extents[num_dims - 1];
    for (int i = 2; i != num_dims; ++i)
      idim[i] = this->desc_.extents[num_dims - i - 1];
  }
  Rf_setAttrib(res, R_DimSymbol, dim);
  Rf_setAttrib(res2, R_DimSymbol, dim);

  if (num_dims <= 2) Rf_copyMatrix(res2, res, TRUE);
  // TODO: else CopyMatrixByRow(res2, res); // for num_dims > 2

  UNPROTECT(3);

  return res2;
}

#endif

//! Matrix<T,0> is a zero-dimensional matrix of type T..
/*!
 * \tparam T value type.
 *
 * Matrix<T,0> is not really a matrix. It stores a single element of
 * type T and can be converted to a reference to that type.
 */
template <typename T>
class Matrix<T, 0> : public MatrixBase<T, 0> {
 public:
  //! @cond Doxygen_Suppress
  using iterator = typename std::array<T, 1>::iterator;
  using const_iterator = typename std::array<T, 1>::const_iterator;
  //! @endcond

  //! construct from an element
  Matrix(const T &x = T{}) : elem_{x} {}
  //! assign from an element
  Matrix &operator=(const T &value) {
    elem_[0] = value;
    return *this;
  }

  //! total number of elements
  std::size_t size() const { return 1; }

  //! "flat" element access
  ///@{
  T *data() { return elem_.data(); }
  const T *data() const { return elem_.data(); }
  ///@}

  //! m() subscripting
  ///@{
  T &operator()() { return elem_[0]; }
  const T &operator()() const { return elem_[0]; }
  ///@}

  //! conversion from Matrix<T, 0> to type T
  ///@{
  operator T &() { return elem_[0]; }
  operator const T &() { return elem_[0]; }
  ///@}

  //! element iterators
  ///@{
  iterator begin() { return elem_.begin(); }
  const_iterator begin() const { return elem_.cbegin(); }
  iterator end() { return elem_.end(); }
  const_iterator end() const { return elem_.end(); }
  ///@}

 private:
  std::array<T, 1> elem_;
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const Matrix<T, 0> &m0) {
  return os << (const T &)m0();
}

_SLAB_END_NAMESPACE

#endif  // _SLAB_MATRIX_MATRIX_H
