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
// matrix_ref.h
// -----------------------------------------------------------------------------
//
#ifndef SLAB_MATRIX_MATRIX_REF_H_
#define SLAB_MATRIX_MATRIX_REF_H_

#include <cstddef>
#include <iterator> // std::forward_iterator_tag
#include "slab/matrix/matrix_base.h"
#include "slab/matrix/matrix.h"

template<typename T, std::size_t N>
class MatrixRefIter;

template<typename T, std::size_t N>
class MatrixRef {
 public:
  static constexpr std::size_t order_ = N;
  using value_type = T;
  using iterator = MatrixRefIter<T, N>;
  using const_iterator = MatrixRefIter<T, N>;

  MatrixRef() = delete;
  MatrixRef(MatrixRef &&) = default;
  MatrixRef &operator=(MatrixRef &&) = default;
  MatrixRef(MatrixRef const &) = default;
  MatrixRef &operator=(MatrixRef const &) = default;
  ~MatrixRef() = default;

  template<typename U>
  MatrixRef &operator=(const Matrix<U, N> &);

  MatrixRef &operator=(MatrixInitializer<T, N>);

  MatrixRef(const MatrixSlice<N> &s, T *p) : desc_{s}, ptr_{p} {}

  template<typename... Args>
  Enable_if<matrix_impl::Requesting_element<Args...>(), T &>
  operator()(Args... args) const;
  template<typename... Args>
  Enable_if<matrix_impl::Requesting_slice<Args...>(), MatrixRef<T, N>>
  operator()(const Args... args) const;

  MatrixRef<T, N - 1> operator[](size_t i) const { return row(i); }
  MatrixRef<T, N - 1> row(size_t n) const;
  MatrixRef<T, N - 1> col(size_t n) const;

  template<typename F>
  MatrixRef &apply(F f);
  MatrixRef &operator=(const T &value);

//    template<template<typename,size_t> class M, typename T2, typename F,
//            typename =Enable_if<Dimensional_Matrix_type<M,T2,N>()>>
//    MatrixRef& apply(const M<T2,N>& m, F f);

  size_t extent(size_t n) const { return desc_.extents[n]; }
  size_t size() const { return desc_.size; }
  size_t rows() const { return desc_.extents[0]; }
  size_t cols() const { return desc_.extents[1]; }
  const MatrixSlice<N> &descriptor() const { return desc_; }

  T *data() const { return ptr_; }

  const_iterator begin() const { return const_iterator::begin(this); }
  const_iterator end() const { return const_iterator::end(this); }

 private:
  MatrixSlice<N> desc_;
  T *ptr_;
};

template<typename T, std::size_t N>
template<typename U>
MatrixRef<T, N> &MatrixRef<T, N>::operator=(const Matrix<U, N> &x) {
  static_assert(Convertible<U, T>(), "MatrixRef =: incompatible element types");
  assert(desc_.extents == x.descriptor().extents);

  std::copy(x.begin(), x.end(), begin());
  return *this;
}

template<typename T, std::size_t N>
MatrixRef<T, N> &MatrixRef<T, N>::operator=(MatrixInitializer<T, N> init) {
  std::array<std::size_t, N> extents = matrix_impl::derive_extents<N>(init);
  assert(extents == desc_.extents);

  auto iter = begin();
  matrix_impl::copy_flat(init, iter);

  return *this;
}

template<typename T, size_t N>
template<typename... Args>
Enable_if<matrix_impl::Requesting_element<Args...>(), T &>
MatrixRef<T, N>::operator()(Args... args) const {
  static_assert(sizeof...(Args) == N,
                "Matrix_ref<T,N>::operator()(size_t...): dimension mismatch");
  assert(matrix_impl::check_bounds(desc_, args...));
  return *(data() + desc_(args...));
}

template<typename T, size_t N>
template<typename... Args>
Enable_if<matrix_impl::Requesting_slice<Args...>(), MatrixRef<T, N>>
MatrixRef<T, N>::operator()(const Args... args) const {
  MatrixSlice<N> d;
  d.start = desc_.start + matrix_impl::do_slice(desc_, d, args...);
  d.size = matrix_impl::compute_size(d.extents);
  return {d, data()};
}

// row

template<typename T, size_t N>
MatrixRef<T, N - 1> MatrixRef<T, N>::row(size_t n) const {
  assert(n < rows());
  MatrixSlice<N - 1> row;
  matrix_impl::slice_dim<0>(n, desc_, row);
  return {row, ptr_};
}

// col

template<typename T, size_t N>
MatrixRef<T, N - 1> MatrixRef<T, N>::col(size_t n) const {
  assert(n < cols());
  MatrixSlice<N - 1> col;
  matrix_impl::slice_dim<1>(n, desc_, col);
  return {col, ptr_};
}

template<typename T>
class MatrixRef<T, 0> {
 public:
  static constexpr size_t order = 0;
  using value_type = T;

  MatrixRef(const MatrixSlice<0> &s, T *p) : ptr_{p + s.start} {}

  T &operator()() { return *ptr_; };
  const T &operator()() const { return *ptr_; }

  operator T &() { return *ptr_; }
  operator const T &() const { return *ptr_; }

 private:
  T *ptr_;
};

template<typename T>
std::ostream &operator<<(std::ostream &os, const MatrixRef<T, 0> &mr0) {
  return os << (const T &) mr0;
}

template<typename T, std::size_t N>
class MatrixRefIter : public std::iterator<std::forward_iterator_tag, T> {
 public:
  static MatrixRefIter begin(const MatrixRef<T, N> *t);
  static MatrixRefIter end(const MatrixRef<T, N> *t);

  MatrixRefIter &operator++();
  T &operator*() { return *(target->data() + target->descriptor().offset(pos)); }
  bool operator!=(const MatrixRefIter &o) { return !(pos == o.pos && ov == o.ov); }

 private:
  MatrixRefIter(const MatrixRef<T, N> *t, bool ov = false)
      : target{t}, ov{ov} { std::fill(pos.begin(), pos.end(), 0); }

  const MatrixRef<T, N> *target;
  std::array<std::size_t, N> pos;
  bool ov;

  template<typename U, size_t NN>
  friend std::ostream &operator<<(std::ostream &os, const MatrixRefIter<U, NN> &iter);
};

template<typename T, size_t N>
MatrixRefIter<T, N> MatrixRefIter<T, N>::begin(const MatrixRef<T, N> *t) {
  MatrixRefIter<T, N> iter{t};
  return iter;
}

template<typename T, size_t N>
MatrixRefIter<T, N> MatrixRefIter<T, N>::end(const MatrixRef<T, N> *t) {
  MatrixRefIter<T, N> iter{t, true};
  return iter;
}

template<typename T, size_t N>
MatrixRefIter<T, N> &MatrixRefIter<T, N>::operator++() {
  for (int i = N - 1; i >= 0; --i) {
    if (++pos[i] < target->descriptor().extents[i]) return *this;
    pos[i] = 0;
  }
  ov = true;
  return *this;
}

template<typename T, size_t N>
std::ostream &operator<<(std::ostream &os, const MatrixRefIter<T, N> &iter) {
  os << "target: " << *iter.target << ", pos: " << iter.pos
     << ", ov: " << iter.ov << std::endl;
  return os;
}

#endif // SLAB_MATRIX_MATRIX_REF_H_
