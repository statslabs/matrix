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

/// @file hermitian_matrix.h
/// @brief Hermitian matrix

#ifndef STATSLABS_MATRIX_HERMITIAN_MATRIX_H_
#define STATSLABS_MATRIX_HERMITIAN_MATRIX_H_

template<typename T, typename TRI>
class HermitianMatrix;

template<typename T, typename TRI>
struct HermitianMatrixElement {
  using value_type = T;

  HermitianMatrixElement(HermitianMatrix<T, TRI> &hm, std::size_t i, std::size_t j, T elem)
      : hm_(hm), i_(i), j_(j), elem_(elem), is_elem_changed_(false) {}

  ~HermitianMatrixElement() { if (is_elem_changed_) hm_.assign_element(i_, j_, elem_); }

  HermitianMatrixElement &operator=(const HermitianMatrixElement &hme) {
    elem_ = hme.elem_;
    is_elem_changed_ = true;
    return *this;
  }

  HermitianMatrixElement &operator=(const T &elem) {
    elem_ = elem;
    is_elem_changed_ = true;
    return *this;
  }

  HermitianMatrixElement &operator+=(const T &elem) {
    elem_ += elem;
    is_elem_changed_ = true;
    return *this;
  }

  HermitianMatrixElement &operator-=(const T &elem) {
    elem_ -= elem;
    is_elem_changed_ = true;
    return *this;
  }

  HermitianMatrixElement &operator*=(const T &elem) {
    elem_ *= elem;
    is_elem_changed_ = true;
    return *this;
  }

  HermitianMatrixElement &operator/=(const T &elem) {
    elem_ /= elem;
    is_elem_changed_ = true;
    return *this;
  }

  bool operator==(const T &elem) const { return elem_ == elem; }
  bool operator!=(const T &elem) const { return elem_ != elem; }
  operator T() const { return elem_; }

 private:
  HermitianMatrix<T, TRI> &hm_;
  std::size_t i_;
  std::size_t j_;
  T elem_;
  bool is_elem_changed_;
};

template<typename T, typename TRI>
class HermitianMatrix {
 public:
  using value_type = T;
  using reference = HermitianMatrixElement<T, TRI>;
  using iterator = typename std::vector<T>::iterator;
  using const_iterator  = typename std::vector<T>::const_iterator;

  HermitianMatrix() = default;
  HermitianMatrix(HermitianMatrix &&) = default;
  HermitianMatrix &operator=(HermitianMatrix &&) = default;
  HermitianMatrix(HermitianMatrix const &) = default;
  HermitianMatrix &operator=(HermitianMatrix const &) = default;

  HermitianMatrix(std::size_t n) : elem_(n * n), desc_(n) {}

  const TRI &descriptor() const { return desc_; }

  //! "flat" element access
  ///@{
  T *data() { return elem_.data(); }
  const T *data() const { return elem_.data(); }
  ///@}

  std::size_t n_rows() const { return desc_.extents[0]; }
  std::size_t n_cols() const { return desc_.extents[1]; }

  reference operator()(std::size_t i, std::size_t j) {
    assert(i < n_rows());
    assert(j < n_cols());

    if (!desc_.other_half(i, j))
      return HermitianMatrixElement<T, TRI>(*this, i, j, *(data() + desc_(i, j)));
    else
      return HermitianMatrixElement<T, TRI>(*this, i, j, std::conj(*(data() + desc_(i, j))));
  }

  const T operator()(std::size_t i, std::size_t j) const {
    assert(i < n_rows());
    assert(j < n_cols());

    if (!desc_.other_half(i, j))
      return *(data() + desc_(i, j));
    else
      return std::conj(*(data() + desc_(j, i)));
  }

  void assign_element(std::size_t i, std::size_t j, const T &val) {
    assert(i < n_rows());
    assert(j < n_cols());

    if (!desc_.other_half(i, j))
      *(data() + desc_(i, j)) = val;
    else
      *(data() + desc_(i, j)) = std::conj(val);
  }

 private:
  std::vector<T> elem_;
  TRI desc_;
};

template<typename T, typename TRI>
std::ostream &operator<<(std::ostream &os, const HermitianMatrix<T, TRI> &m) {
  for (std::size_t i = 0; i != m.n_rows(); ++i) {
    for (std::size_t j = 0; j != m.n_cols(); ++j) {
      os << T(m(i, j)) << "\t";
    }
    os << std::endl;
  }

  return os << std::endl;
}

#endif // STATSLABS_MATRIX_HERMITIAN_MATRIX_H_
