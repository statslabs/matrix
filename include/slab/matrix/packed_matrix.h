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

/// @file packed_matrix.h
/// @brief Packed matrix

#ifndef SLAB_MATRIX_PACKED_MATRIX_H_
#define SLAB_MATRIX_PACKED_MATRIX_H_

#include <cstddef>

#include <array>
#include <complex>
#include <iostream>
#include <vector>

#include "slab/matrix/error.h"
#include "slab/matrix/traits.h"

namespace slab {

// Triangular matrix type
struct lower_tag {};
struct upper_tag {};
struct unit_lower_tag : public lower_tag {};
struct unit_upper_tag : public upper_tag {};

struct upper {
  using triangular_type = upper_tag;

  upper();
  upper(std::size_t n) {
    extents[0] = n;
    extents[1] = n;
    size = std::size_t((1 + n) * n / 2);
  }

  inline bool other_half(std::size_t i, std::size_t j) const { return i > j; }

  inline std::size_t operator()(std::size_t i, std::size_t j) const {
    return (2 * extents[0] - i + 1) * i / 2 + (j - i);
  }

  std::size_t size;
  std::array<std::size_t, 2> extents;
};

struct unit_upper : public upper {
  using triangular_type = unit_upper_tag;
  unit_upper(std::size_t n) : upper(n) { size -= n; }
};

struct lower {
  using triangular_type = lower_tag;

  lower();
  lower(std::size_t n) {
    extents[0] = n;
    extents[1] = n;
    size = std::size_t((1 + n) * n / 2);
  }

  inline bool other_half(std::size_t i, std::size_t j) const { return j > i; }

  inline std::size_t operator()(std::size_t i, std::size_t j) const {
    return (1 + i) * i / 2 + j;
  }

  std::size_t size;
  std::array<std::size_t, 2> extents;
};

struct unit_lower : public lower {
  using triangular_type = unit_lower_tag;
  unit_lower(std::size_t n) : lower(n) { size -= n; }
};

template <typename T>
struct is_upper : public std::false_type {};

template <>
struct is_upper<upper> : public std::true_type {};

template <typename T>
struct is_lower : public std::false_type {};

template <>
struct is_lower<lower> : public std::true_type {};

template <typename T, typename TRI>
class PackedMatrix;

template <typename T, typename TRI>
struct PackedMatrixElement {
  using value_type = T;

  PackedMatrixElement(PackedMatrix<T, TRI> &hm, std::size_t i, std::size_t j,
                      T elem)
      : hm_(hm), i_(i), j_(j), elem_(elem), is_elem_changed_(false) {}

  ~PackedMatrixElement() {
    if (is_elem_changed_) {
      hm_.assign_element(i_, j_, elem_);
    }
  }

  PackedMatrixElement &operator=(const PackedMatrixElement &hme) {
    elem_ = hme.elem_;
    is_elem_changed_ = true;
    return *this;
  }

  PackedMatrixElement &operator=(const T &elem) {
    elem_ = elem;
    is_elem_changed_ = true;
    return *this;
  }

  PackedMatrixElement &operator+=(const T &elem) {
    elem_ += elem;
    is_elem_changed_ = true;
    return *this;
  }

  PackedMatrixElement &operator-=(const T &elem) {
    elem_ -= elem;
    is_elem_changed_ = true;
    return *this;
  }

  PackedMatrixElement &operator*=(const T &elem) {
    elem_ *= elem;
    is_elem_changed_ = true;
    return *this;
  }

  PackedMatrixElement &operator/=(const T &elem) {
    elem_ /= elem;
    is_elem_changed_ = true;
    return *this;
  }

  bool operator==(const T &elem) const { return elem_ == elem; }
  bool operator!=(const T &elem) const { return elem_ != elem; }
  operator T() const { return elem_; }

 private:
  PackedMatrix<T, TRI> &hm_;
  std::size_t i_;
  std::size_t j_;
  T elem_;
  bool is_elem_changed_;
};

template <typename T, typename TRI>
class PackedMatrix {
 public:
  using value_type = T;
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;
  using triangular_type = typename TRI::triangular_type;
  using reference = PackedMatrixElement<T, TRI>;

  PackedMatrix() = default;
  PackedMatrix(PackedMatrix &&) = default;
  PackedMatrix &operator=(PackedMatrix &&) = default;
  PackedMatrix(PackedMatrix const &) = default;
  PackedMatrix &operator=(PackedMatrix const &) = default;

  PackedMatrix(std::size_t n) : elems_(n * n), desc_(n) {}
  PackedMatrix(std::size_t n, const std::vector<T> &v) : desc_(n) {
    if (v.size() != desc_.size)
      err_quit(
          "Fail to construct a packed matrix, the size of vector provided is "
          "incorrect");

    init(v, triangular_type());
  }
  PackedMatrix(std::size_t n, const std::initializer_list<T> &il) : desc_(n) {
    if (il.size() != desc_.size)
      err_quit(
          "Fail to construct a packed matrix, the size of vector provided is "
          "incorrect");

    std::vector<T> v(il);
    init(v, triangular_type());
  }

  const TRI &descriptor() const { return desc_; }

  //! "flat" element access
  ///@{
  T *data() { return elems_.data(); }
  const T *data() const { return elems_.data(); }
  ///@}

  std::size_t n_rows() const { return desc_.extents[0]; }
  std::size_t n_cols() const { return desc_.extents[1]; }

  reference operator()(std::size_t i, std::size_t j) {
    assert(i < n_rows());
    assert(j < n_cols());

    if (!this->desc_.other_half(i, j))
      return PackedMatrixElement<T, TRI>(*this, i, j,
                                         *(this->data() + this->desc_(i, j)));
    else
      return PackedMatrixElement<T, TRI>(
          *this, i, j,
          element_in_other_half(*(this->data() + this->desc_(j, i))));
  }

  const T operator()(std::size_t i, std::size_t j) const {
    assert(i < n_rows());
    assert(j < n_cols());

    if (!this->desc_.other_half(i, j))
      return *(this->data() + this->desc_(i, j));
    else
      return element_in_other_half(*(this->data() + this->desc_(j, i)));
  }

  void assign_element(std::size_t i, std::size_t j, const T &val) {
    assert(i < n_rows());
    assert(j < n_cols());

    if (!this->desc_.other_half(i, j))
      *(this->data() + this->desc_(i, j)) = val;
    else
      update_element_in_base_half(*(this->data() + this->desc_(j, i)), val);
  }

 protected:
  std::vector<T> elems_;
  TRI desc_;

  // Given the element in the base half,
  // return the value of element in the corresponding other half
  virtual T element_in_other_half(const T &val) const = 0;

  // Given the value of element in the other half,
  // update the value of element in the corresponding base half
  virtual void update_element_in_base_half(T &elem, const T &val) = 0;

 private:
  void init(const std::vector<T> &v, upper_tag) {
    elems_.assign(v.begin(), v.end());
  }
  void init(const std::vector<T> &v, lower_tag) {
    elems_.assign(v.begin(), v.end());
  }
  void init(const std::vector<T> &v, unit_upper_tag) {
    auto iter = v.begin();
    for (std::size_t j = 0; j != desc_.extents[1]; ++j) {
      for (std::size_t i = j; i != desc_.extents[0]; ++i) {
        if (i == j)
          elems_.push_back(T{1});
        else
          elems_.push_back(*iter++);
      }
    }
  }
  void init(const std::vector<T> &v, unit_lower_tag) {
    auto iter = v.begin();
    for (std::size_t i = 0; i != desc_.extents[0]; ++i) {
      for (std::size_t j = 0; j <= i; ++j) {
        if (i == j)
          elems_.push_back(T{1});
        else
          elems_.push_back(*iter++);
      }
    }
  }
};

template <typename T, typename TRI>
class SymmetricMatrix : public PackedMatrix<T, TRI> {
 public:
  SymmetricMatrix(std::size_t n) : PackedMatrix<T, TRI>{n} {}
  SymmetricMatrix(std::size_t n, const std::vector<T> &v)
      : PackedMatrix<T, TRI>(n, v) {}

 private:
  T element_in_other_half(const T &val) const override { return val; }
  void update_element_in_base_half(T &elem, const T &val) override {
    elem = val;
  }
};

template <typename T, typename TRI>
class TriangularMatrix : public PackedMatrix<T, TRI> {
 public:
  TriangularMatrix(std::size_t n) : PackedMatrix<T, TRI>{n} {}
  TriangularMatrix(std::size_t n, const std::vector<T> &v)
      : PackedMatrix<T, TRI>(n, v) {}

 private:
  T element_in_other_half(const T &val) const override {
    ignore(val);
    return 0;
  }
  void update_element_in_base_half(T &elem, const T &val) override {
    ignore(elem);
    ignore(val);
  }
};

template <typename T, typename TRI>
class HermitianMatrix : public PackedMatrix<T, TRI> {
 public:
  HermitianMatrix(std::size_t n) : PackedMatrix<T, TRI>{n} {}

 private:
  T element_in_other_half(const T &val) const override {
    return std::conj(val);
  }
  void update_element_in_base_half(T &elem, const T &val) override {
    elem = std::conj(val);
  }
};

template <typename T, typename TRI>
std::ostream &operator<<(std::ostream &os, const PackedMatrix<T, TRI> &m) {
  for (std::size_t i = 0; i != m.n_rows(); ++i) {
    for (std::size_t j = 0; j != m.n_cols(); ++j) {
      os << m(i, j) << "\t";
    }
    os << std::endl;
  }

  return os << std::endl;
}

}  // namespace slab

#endif  // SLAB_MATRIX_PACKED_MATRIX_H
