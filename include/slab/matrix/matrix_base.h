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
// matrix_base.h
// -----------------------------------------------------------------------------
//
#ifndef SLAB_MATRIX_MATRIX_BASE_H_
#define SLAB_MATRIX_MATRIX_BASE_H_

template<typename T, std::size_t N>
class MatrixBase {
 public:
  static constexpr std::size_t order_ = N;
  using value_type = T;

  MatrixBase() = default;
  MatrixBase(MatrixBase &&) = default;
  MatrixBase &operator=(MatrixBase &&) = default;
  MatrixBase(MatrixBase const &) = default;
  MatrixBase &operator=(MatrixBase const &) = default;

  template<typename... Exts>
  explicit MatrixBase(Exts... exts) : desc_{exts...} {}

  explicit MatrixBase(const MatrixSlice<N> &ms) : desc_{ms} {}

  // number of dimensions
  static constexpr std::size_t order() { return order_; }
  // #elements in the nth dimension
  std::size_t extent(std::size_t n) const {
    assert(n < order_);
    return desc_.extents[n];
  }
  // total number of elements
  virtual std::size_t size() const = 0;
  // the slice defining subscripting
  const MatrixSlice<N> &descriptor() const { return desc_; }

  virtual T *data() = 0;
  virtual const T *data() const = 0;

  std::size_t n_rows() const { return desc_.extents[0]; }
  std::size_t n_cols() const { return desc_.extents[1]; }

  // m(i,j,k) subscripting with integers
  template<typename... Args>
  T &operator()(Args... args);

  template<typename... Args>
  const T &operator()(Args... args) const;

 protected:
  MatrixSlice<N> desc_;
};

template<typename T, std::size_t N>
template<typename... Args>
T &MatrixBase<T, N>::operator()(Args... args) {
  assert(matrix_impl::check_bounds(this->desc_, args...));
  return *(data() + this->desc_(args...));
}

template<typename T, std::size_t N>
template<typename... Args>
const T &MatrixBase<T, N>::operator()(Args... args) const {
  assert(matrix_impl::check_bounds(this->desc_, args...));
  return *(data() + this->desc_(args...));
}

#endif // SLAB_MATRIX_MATRIX_BASE_H_
