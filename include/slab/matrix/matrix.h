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

/// @file matrix.h
/// @brief A Matrix template

#ifndef SLAB_MATRIX_MATRIX_H_
#define SLAB_MATRIX_MATRIX_H_

#include <cassert>
#include <cstddef>

#include <fstream>
#include <string>
#include <vector>

#include "slab/matrix/matrix_base.h"
#include "slab/matrix/matrix_ref.h"
#include "slab/matrix/matrix_slice.h"
#include "slab/matrix/packed_matrix.h"
#include "slab/matrix/support.h"

namespace slab {

template <typename T>
Matrix<T, 2> transpose(const MatrixBase<T, 1> &a);

template <typename T>
Matrix<T, 2> transpose(const MatrixBase<T, 2> &a);

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
  template <typename U, std::size_t NN = N, typename = Enable_if<(NN > 1)>>
  Matrix(std::initializer_list<U>) = delete;
  template <typename U, std::size_t NN = N, typename = Enable_if<(NN > 1)>>
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
  MatrixRef<T, N - 1> col(size_t n);
  MatrixRef<const T, N - 1> col(size_t n) const;
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
  //! @cond Doxygen_Suppress
  template <typename F>
  Matrix &apply(F f);  // f(x) for every element x

  // f(x, mx) for corresponding elements of *this and m
  template <typename M, typename F>
  Enable_if<Matrix_type<M>(), Matrix &> apply(const M &m, F f);

  Matrix &operator=(const T &value);   // assignment with scalar

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

  Matrix operator-() const;
  //! @endcond

  // ---------------------------------------------
  // Specializations for Matrix<T, 1>/Matrix<T, 2>
  // ---------------------------------------------

 public:
  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 1)>>
  Matrix(const Matrix<U, 2> &x)
      : MatrixBase<T, N>{x.n_rows()}, elems_{x.begin(), x.end()} {
    static_assert(Convertible<U, T>(),
                  "Matrix constructor: incompatible element types");
    assert(x.n_cols() == 1);
  }

  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 1)>>
  Matrix(const MatrixRef<U, 2> &x)
      : MatrixBase<T, N>{x.n_rows()}, elems_{x.begin(), x.end()} {
    static_assert(Convertible<U, T>(),
                  "Matrix constructor: incompatible element types");
    assert(x.n_cols() == 1);
  }

  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 1)>>
  Matrix &operator=(const Matrix<U, 2> &x) {
    static_assert(Convertible<U, T>(), "Matrix =: incompatible element types");
    assert(x.n_cols() == 1);

    this->desc_.size = x.descriptor().size;
    this->desc_.start = 0;
    this->desc_.extents[0] = x.n_rows();
    this->desc_.strides[0] = 1;

    elems_.assign(x.begin(), x.end());

    return *this;
  }

  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 1)>>
  Matrix &operator=(const MatrixRef<U, 2> &x) {
    static_assert(Convertible<U, T>(), "Matrix =: incompatible element types");
    assert(x.n_cols() == 1);

    this->desc_.size = x.descriptor().size;
    this->desc_.start = 0;
    this->desc_.extents[0] = x.n_rows();
    this->desc_.strides[0] = 1;

    elems_.assign(x.begin(), x.end());

    return *this;
  }

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

  template <std::size_t NN = N, typename = Enable_if<(NN == 1)>>
  std::vector<T> std_vec() const {
    return std::vector<T>(begin(), end());
  }

 public:
  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 2)>>
  Matrix(const Matrix<U, 1> &x)
      : MatrixBase<T, N>{x.n_rows(), 1}, elems_{x.begin(), x.end()} {
    static_assert(Convertible<U, T>(),
                  "Matrix constructor: incompatible element types");
  }

  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 2)>>
  Matrix(const MatrixRef<U, 1> &x)
      : MatrixBase<T, N>{x.n_rows(), 1}, elems_{x.begin(), x.end()} {
    static_assert(Convertible<U, T>(),
                  "Matrix constructor: incompatible element types");
  }

  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 2)>>
  Matrix &operator=(const Matrix<U, 1> &x) {
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

  template <typename U, std::size_t NN = N, typename = Enable_if<(NN == 2)>>
  Matrix &operator=(const MatrixRef<U, 1> &x) {
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

  template <typename U, typename TRI, std::size_t NN = N,
            typename = Enable_if<(NN == 2)>>
  Matrix(const SymmetricMatrix<U, TRI> &x)
      : MatrixBase<T, N>{x.n_rows(), x.n_cols()} {
    static_assert(Convertible<U, T>(),
                  "Matrix constructor: incompatible element types");
    for (std::size_t i = 0; i != x.n_rows(); ++i) {
      for (std::size_t j = 0; j != x.n_cols(); ++j) {
        elems_.push_back(x(i, j));
      }
    }
  }

  template <typename U, typename TRI, std::size_t NN = N,
            typename = Enable_if<(NN == 2)>>
  Matrix(const TriangularMatrix<U, TRI> &x)
      : MatrixBase<T, N>{x.n_rows(), x.n_cols()} {
    static_assert(Convertible<U, T>(),
                  "Matrix constructor: incompatible element types");
    for (std::size_t i = 0; i != x.n_rows(); ++i) {
      for (std::size_t j = 0; j != x.n_cols(); ++j) {
        elems_.push_back(x(i, j));
      }
    }
  }

  template <typename U, typename TRI, std::size_t NN = N,
            typename = Enable_if<(NN == 2)>>
  Matrix(const HermitianMatrix<U, TRI> &x)
      : MatrixBase<T, N>{x.n_rows(), x.n_cols()} {
    static_assert(Convertible<U, T>(),
                  "Matrix constructor: incompatible element types");
    for (std::size_t i = 0; i != x.n_rows(); ++i) {
      for (std::size_t j = 0; j != x.n_cols(); ++j) {
        elems_.push_back(x(i, j));
      }
    }
  }

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

  template <std::size_t NN = N, typename = Enable_if<(NN == 1) || (NN == 2)>>
  Matrix<T, 2> t() const {
    return transpose(*this);
  }

  // ----------------------------------------
  // More member function for matrix template
  // ----------------------------------------

 public:
  void clear();

  bool empty() const { return begin() == end(); }
  bool is_empty() const { return empty(); }

  //  void save(const std::string &filename) {
  //    std::ostream os(filename);
  //  }
  void load(const std::string &filename) {
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
      this->desc_.size = matrix_impl::compute_strides(this->desc_.extents,
                                                      this->desc_.strides);

      std::istream_iterator<T> in(is), end;
      elems_.assign(in, end);
    } else {
      std::cout << "Fail to open the file" << std::endl;
    }
  }
};

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
Matrix<T, N> Matrix<T, N>::operator-() const {
  Matrix<T, N> res = *this;
  return res.apply([&](T &a) { a = -a; });
}

template <typename T, std::size_t N>
void Matrix<T, N>::clear() {
  this->desc_.clear();
  elems_.clear();
}

template <typename T>
class Matrix<T, 0> : public MatrixBase<T, 0> {
 public:
  using iterator = typename std::array<T, 1>::iterator;
  using const_iterator = typename std::array<T, 1>::const_iterator;

  Matrix(const T &x = T{}) : elem_{x} {}

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

  T &operator()() { return elem_[0]; }
  const T &operator()() const { return elem_[0]; }

  operator T &() { return elem_[0]; }
  operator const T &() { return elem_[0]; }

  iterator begin() { return elem_.begin(); }
  const_iterator begin() const { return elem_.cbegin(); }
  iterator end() { return elem_.end(); }
  const_iterator end() const { return elem_.end(); }

 private:
  std::array<T, 1> elem_;
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const Matrix<T, 0> &m0) {
  return os << (const T &)m0();
}

template <typename M, typename... Args>
Enable_if<Matrix_type<M>(), M> zeros(Args... args) {
  assert(M::order() == sizeof...(args));
  M res(args...);
  res = 0;

  return res;
}

template <typename M, typename... Args>
Enable_if<Matrix_type<M>(), M> ones(Args... args) {
  assert(M::order() == sizeof...(args));
  M res(args...);
  res = 1;

  return res;
}

template <typename M, typename... Args>
Enable_if<Matrix_type<M>(), M> eye(std::size_t i, std::size_t j) {
  assert(M::order() == 2);
  M res(i, j);
  res.diag() = 1;

  return res;
}

template <typename T>
Matrix<T, 2> transpose(const MatrixBase<T, 1> &a) {
  Matrix<T, 2> res(1, a.n_rows());
  for (std::size_t i = 0; i < a.n_rows(); ++i) {
    res(0, i) = a(i);
  }

  return res;
}

template <typename T>
Matrix<T, 2> transpose(const MatrixBase<T, 2> &a) {
  Matrix<T, 2> res(a.n_cols(), a.n_rows());
  for (std::size_t i = 0; i < a.n_rows(); ++i) {
    for (std::size_t j = 0; j < a.n_cols(); ++j) {
      res(j, i) = a(i, j);
    }
  }

  return res;
}

}  // namespace slab

#endif  // SLAB_MATRIX_MATRIX_H_
