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
// matrix_slice.h
// -----------------------------------------------------------------------------
//
#ifndef SLAB_MATRIX_MATRIX_SLICE_H_
#define SLAB_MATRIX_MATRIX_SLICE_H_

#include <cstddef>
#include <algorithm>
#include <array>
#include <initializer_list>
#include "slab/matrix/traits.h"

// A matrix slice specifies the N-dimensional matrix properties of a contiguous
// region of memory. The slice is primarily described by 3 parameters:
// A sequence of extents, a sequence of strides, and a starting offset. An
// extent describes an array bound in a particular dimension. A stride defines
// the distance between sub-matrices in a particular dimension. This class
// also maintains the total number of elements in the slice, which is just the
// product of dimensions.
template<std::size_t N>
struct MatrixSlice {

  MatrixSlice();                               // an empty matrix: no elements

  MatrixSlice(std::size_t s, std::initializer_list<std::size_t> exts); // extents
  MatrixSlice(std::size_t s, std::initializer_list<std::size_t> exts,  // extents and strides
              std::initializer_list<std::size_t> strs);
  MatrixSlice(const std::array<std::size_t, N> &exts);

  template<typename... Dims>
  MatrixSlice(Dims... dims);                   // N extents

  template<typename... Dims>
  std::size_t operator()(Dims... dims) const;  // calculate index from a set of subscripts

  std::size_t offset(const std::array<std::size_t, N> &pos) const;

  std::size_t size;                   // total number of elements
  std::size_t start;                  // starting offset
  std::array<std::size_t, N> extents; // number of elements in each dimension
  std::array<std::size_t, N> strides; // offsets between elements in each dimension
};

template<size_t N>
MatrixSlice<N>::MatrixSlice()
    :size{0}, start{0} {
  std::fill(extents.begin(), extents.end(), 0);
  std::fill(strides.begin(), strides.end(), 1);
}

template<std::size_t N>
MatrixSlice<N>::MatrixSlice(std::size_t s,
                            std::initializer_list<std::size_t> exts)
    : start(s) {
  assert(exts.size() == N);
  std::copy(exts.begin(), exts.end(), extents.begin());
  size = matrix_impl::compute_strides(extents, strides);
}

template<std::size_t N>
MatrixSlice<N>::MatrixSlice(std::size_t s,
                            std::initializer_list<std::size_t> exts,
                            std::initializer_list<std::size_t> strs)
    : start(s) {
  assert(exts.size() == N);
  std::copy(exts.begin(), exts.end(), extents.begin());
  std::copy(strs.begin(), strs.end(), strides.begin());
  size = matrix_impl::compute_size(extents);
}

template<std::size_t N>
MatrixSlice<N>::MatrixSlice(const std::array<std::size_t, N> &exts)
    : start{0}, extents{exts} {
  assert(exts.size() == N);
  size = matrix_impl::compute_strides(extents, strides);
}

template<std::size_t N>
template<typename... Dims>
MatrixSlice<N>::MatrixSlice(Dims... dims)
    : start{0} {
  static_assert(sizeof...(Dims) == N,
                "Matrix_slice<N>::Matrix_slice(Dims...): dimension mismatch");

  std::size_t args[N]{std::size_t(dims)...};
  std::copy(std::begin(args), std::end(args), extents.begin());
  size = matrix_impl::compute_strides(extents, strides);
}

template<std::size_t N>
template<typename... Dims>
std::size_t MatrixSlice<N>::operator()(Dims... dims) const {
  static_assert(sizeof...(Dims) == N, "Matrix_slice<N>::operator(): dimension mismatch");
  std::size_t args[N]{std::size_t(dims)...};  // copy arguments into an array
  return start + std::inner_product(args, args + N, strides.begin(), std::size_t{0});
}

template<size_t N>
std::size_t MatrixSlice<N>::offset(const std::array<std::size_t, N> &pos) const {
  assert(pos.size() == N);
  return start + std::inner_product(pos.begin(), pos.end(), strides.begin(), size_t{0});
}

template<size_t N>
bool same_extents(const MatrixSlice<N>& a, const MatrixSlice<N>& b)
{
  return a.extents == b.extents;
}


template<std::size_t N>
std::ostream &operator<<(std::ostream &os, const std::array<std::size_t, N> &a) {
  for (auto x : a) os << x << ' ';
  return os;
}

template<std::size_t N>
std::ostream &operator<<(std::ostream &os, const MatrixSlice<N> &ms) {
  os << "size: " << ms.size
     << ", start: " << ms.start
     << ", extents: " << ms.extents
     << ", strides: " << ms.strides;
  return os;
}

#endif // SLAB_MATRIX_MATRIX_SLICE_H_