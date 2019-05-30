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

/// @file matrix_base.h
/// @brief Matrix base template

#ifndef _SLAB_MATRIX_MATRIX_BASE_H
#define _SLAB_MATRIX_MATRIX_BASE_H

#include <cassert>
#include <cstddef>
#include <cstdio>

#include <iostream>

#include "slab/__config"

#include "slab/matrix/matrix_slice.h"
#include "slab/matrix/support.h"
#include "slab/matrix/traits.h"

_SLAB_BEGIN_NAMESPACE

template <typename T, std::size_t N>
class MatrixBase {
 public:
  static constexpr std::size_t order_ = N;  // number of dimensions
  using value_type = T;

  MatrixBase() = default;
  MatrixBase(MatrixBase &&) = default;  // move
  MatrixBase &operator=(MatrixBase &&) = default;
  MatrixBase(const MatrixBase &) = default;  // copy
  MatrixBase &operator=(const MatrixBase &) = default;
  ~MatrixBase() = default;

  //! specify the extents
  template <typename... Exts>
  explicit MatrixBase(Exts... exts) : desc_{exts...} {}

  explicit MatrixBase(const MatrixSlice<N> &ms) : desc_{ms} {}

  //! number of dimensions
  static constexpr std::size_t order() { return order_; }
  //! #elements in the nth dimension
  std::size_t extent(std::size_t n) const {
    assert(n < order_);
    return desc_.extents[n];
  }
  //! total number of elements
  virtual std::size_t size() const = 0;
  //! the slice defining subscripting
  const MatrixSlice<N> &descriptor() const { return desc_; }

  //! "flat" element access
  ///@{
  virtual T *data() = 0;
  virtual const T *data() const = 0;
  ///@}

  std::size_t n_rows() const { return desc_.extents[0]; }
  std::size_t n_cols() const { return desc_.extents[1]; }

  // m(i,j,k) subscripting with integers
  template <typename... Args>
  Enable_if<matrix_impl::Requesting_element<Args...>(), T &> operator()(
      Args... args);

  template <typename... Args>
  Enable_if<matrix_impl::Requesting_element<Args...>(), const T &> operator()(
      Args... args) const;

 protected:
  MatrixSlice<N> desc_;  // slice defining extents in the N dimensions
};

template <typename T, std::size_t N>
template <typename... Args>
Enable_if<matrix_impl::Requesting_element<Args...>(), T &> MatrixBase<T, N>::
operator()(Args... args) {
  assert(matrix_impl::check_bounds(this->desc_, args...));
  return *(data() + this->desc_(args...));
}

template <typename T, std::size_t N>
template <typename... Args>
Enable_if<matrix_impl::Requesting_element<Args...>(), const T &>
MatrixBase<T, N>::operator()(Args... args) const {
  assert(matrix_impl::check_bounds(this->desc_, args...));
  return *(data() + this->desc_(args...));
}

template <typename M>
Enable_if<Matrix_type<M>(), std::ostream &> operator<<(std::ostream &os,
                                                       const M &m) {
  os << '{';
  for (std::size_t i = 0; i != m.n_rows(); ++i) {
    os << m[i];
    if (i + 1 != m.n_rows()) os << ',';
  }
  return os << '}';
}

template <typename M>
void raw_print(const M &m) {
  for (auto iter = m.begin(); iter != m.end(); ++iter) {
    printf(" %6.2f", *iter);
  }
}

_SLAB_END_NAMESPACE

#endif  // _SLAB_MATRIX_MATRIX_BASE_H
