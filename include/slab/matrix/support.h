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
//
// -----------------------------------------------------------------------------
// support.h
// -----------------------------------------------------------------------------
//
#ifndef SLAB_MATRIX_SUPPORT_H_
#define SLAB_MATRIX_SUPPORT_H_

#include <cstddef>

#include <algorithm>
#include <initializer_list>

#include "slab/matrix/matrix_fwd.h"
#include "slab/matrix/matrix_slice.h"
#include "slab/matrix/traits.h"

namespace slab {
namespace matrix_impl {

template <std::size_t N, typename List>
bool check_non_jagged(const List &list);

// Describes the structure of a nested std::initializer_list and
// has MatrixInit<T, N - 1> as its member type.
template <typename T, std::size_t N>
struct MatrixInit {
  using type = std::initializer_list<typename MatrixInit<T, N - 1>::type>;
};

// The N == 1 is special. That's where we got to the (most deeply nested).
// std::initialize_list<T>
template <typename T>
struct MatrixInit<T, 1> {
  using type = std::initializer_list<T>;
};

// To avoid surprices, we define N = 0 to be an error
template <typename T>
struct MatrixInit<T, 0>;  // undefined on purpose

template <std::size_t N, typename I, typename List>
Enable_if<(N == 1), void> add_extents(I &first, const List &list) {
  *first = list.size();
}

template <std::size_t N, typename I, typename List>
Enable_if<(N > 1), void> add_extents(I &first, const List &list) {
  assert(check_non_jagged<N>(list));
  *first++ = list.size();
  add_extents<N - 1>(first, *list.begin());
}

// Determines the shape of the Matrix:
//   + Checks that tree really is N deep;
//   + Checks that each row (sub-initializer_list) has the same number of
//   elements;
//   + Sets the extents of each row.
template <std::size_t N, typename List>
std::array<std::size_t, N> derive_extents(const List &list) {
  std::array<std::size_t, N> a;
  auto f = a.begin();
  add_extents<N>(f, list);
  return a;
}

// Checks that all rows have the same number of elements.
template <std::size_t N, typename List>
bool check_non_jagged(const List &list) {
  auto i = list.begin();
  for (auto j = i + 1; j != list.end(); ++j) {
    if (derive_extents<N - 1>(*i) != derive_extents<N - 1>(*j)) return false;
  }
  return true;
}

template <std::size_t N>
std::size_t compute_strides(const std::array<std::size_t, N> &exts,
                            std::array<std::size_t, N> &strs) {
  std::size_t st = 1;
  for (int i = N - 1; i >= 0; --i) {
    strs[i] = st;
    st *= exts[i];
  }
  return st;
}

template <std::size_t N>
std::size_t compute_size(const std::array<size_t, N> &exts) {
  return std::accumulate(exts.begin(), exts.end(), 1,
                         std::multiplies<std::size_t>{});
}

template <typename T, typename Vec>
void add_list(const T *first, const T *last, Vec &vec) {
  vec.insert(vec.end(), first, last);
}

template <typename T, typename Vec>
// nested initializer_lists
void add_list(const std::initializer_list<T> *first,
              const std::initializer_list<T> *last, Vec &vec) {
  for (; first != last; ++first) add_list(first->begin(), first->end(), vec);
}

// Copies the elements of the tree of std::initializer_list s for a Matrix<T,
// N>.
template <typename T, typename Vec>
void insert_flat(std::initializer_list<T> list, Vec &vec) {
  add_list(list.begin(), list.end(), vec);
}

template <typename T, typename Iter>
void copy_list(const T *first, const T *last, Iter &iter) {
  iter = std::copy(first, last, iter);
}

template <typename T, typename Iter>
void copy_list(const std::initializer_list<T> *first,
               const std::initializer_list<T> *last, Iter &it) {
  for (; first != last; ++first) copy_list(first->begin(), first->end(), it);
}

template <typename T, typename Iter>
void copy_flat(std::initializer_list<T> list, Iter &iter) {
  copy_list(list.begin(), list.end(), iter);
}

template <std::size_t I, std::size_t N>
void slice_dim(std::size_t offset, const MatrixSlice<N> &desc,
               MatrixSlice<N - 1> &row) {
  row.start = desc.start;

  int j = (int)N - 2;
  for (int i = N - 1; i >= 0; --i) {
    if (i == I)
      row.start += desc.strides[i] * offset;
    else {
      row.extents[j] = desc.extents[i];
      row.strides[j] = desc.strides[i];
      --j;
    }
  }

  row.size = compute_size(row.extents);
}

template <typename... Args>
constexpr bool Requesting_element() {
  return All(Convertible<Args, size_t>()...);
}

template <typename... Args>
constexpr bool Requesting_slice() {
  return All((Convertible<Args, size_t>() || Same<Args, slice>())...) &&
         Some(Same<Args, slice>()...);
}

template <std::size_t N, typename... Dims>
bool check_bounds(const MatrixSlice<N> &ms, Dims... dims) {
  std::size_t indexes[N]{std::size_t(dims)...};
  return std::equal(indexes, indexes + N, ms.extents.begin(),
                    std::less<std::size_t>{});
}

template <std::size_t NRest, std::size_t N>
std::size_t do_slice_dim(const MatrixSlice<N> &os, MatrixSlice<N> &ns,
                         std::size_t s) {
  std::size_t i = N - NRest;
  ns.strides[i] = os.strides[i];
  ns.extents[i] = 1;
  return s * ns.strides[i];
}

template <std::size_t NRest, std::size_t N>
std::size_t do_slice_dim(const MatrixSlice<N> &os, MatrixSlice<N> &ns,
                         slice s) {
  std::size_t i = N - NRest;
  ns.strides[i] = s.stride * os.strides[i];
  ns.extents[i] = (s.length == size_t(-1))
                      ? (os.extents[i] - s.start + s.stride - 1) / s.stride
                      : s.length;
  return s.start * os.strides[i];
}

template <std::size_t N>
std::size_t do_slice_dim2(const MatrixSlice<N> &os, MatrixSlice<N> &ns, slice s,
                          std::size_t NRest) {
  std::size_t i = N - NRest;
  ns.strides[i] = s.stride * os.strides[i];
  ns.extents[i] = (s.length == size_t(-1))
                      ? (os.extents[i] - s.start + s.stride - 1) / s.stride
                      : s.length;
  return s.start * os.strides[i];
}

template <std::size_t N>
std::size_t do_slice(const MatrixSlice<N> &os, MatrixSlice<N> &ns) {
  return 0;
}

template <std::size_t N, typename T, typename... Args>
std::size_t do_slice(const MatrixSlice<N> &os, MatrixSlice<N> &ns, const T &s,
                     const Args &... args) {
  std::size_t m = do_slice_dim<sizeof...(Args) + 1>(os, ns, s);
  std::size_t n = do_slice(os, ns, args...);
  return m + n;
}

}  // namespace matrix_impl

template <typename T, std::size_t N>
using MatrixInitializer = typename matrix_impl::MatrixInit<T, N>::type;

}  // namespace slab

#endif  // SLAB_MATRIX_SUPPORT_H_
