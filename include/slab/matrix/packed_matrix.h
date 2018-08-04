#ifndef SLAB_MATRIX_PACKED_MATRIX_H_
#define SLAB_MATRIX_PACKED_MATRIX_H_

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

  inline bool other_half(std::size_t i, std::size_t j) const {
    return i > j;
  }

  inline std::size_t operator()(std::size_t i, std::size_t j) const {
    return (2 * extents[0] - i + 1) * i / 2 + (j - i);
  }

  std::size_t size;
  std::array<std::size_t, 2> extents;
};

struct lower {
  using triangular_type = lower_tag;

  lower();
  lower(std::size_t n) {
    extents[0] = n;
    extents[1] = n;
    size = std::size_t((1 + n) * n / 2);
  }

  inline bool other_half(std::size_t i, std::size_t j) const {
    return j > i;
  }

  inline std::size_t operator()(std::size_t i, std::size_t j) const {
    return (1 + i) * i / 2 + j;
  }

  std::size_t size;
  std::array<std::size_t, 2> extents;
};

template<typename T>
struct is_upper : public std::false_type {};

template<>
struct is_upper<upper> : public std::true_type {};

template<typename T>
struct is_lower : public std::false_type {};

template<>
struct is_lower<lower> : public std::true_type {};

template<typename T, typename TRI>
class PackedMatrix {
 public:
  using value_type = T;
  using iterator = typename std::vector<T>::iterator;
  using const_iterator  = typename std::vector<T>::const_iterator;

  PackedMatrix() = default;
  PackedMatrix(PackedMatrix &&) = default;
  PackedMatrix &operator=(PackedMatrix &&) = default;
  PackedMatrix(PackedMatrix const &) = default;
  PackedMatrix &operator=(PackedMatrix const &) = default;

  PackedMatrix(std::size_t n) : elem_(n * n), desc_(n) {}

  const TRI &descriptor() const { return desc_; }

  //! "flat" element access
  ///@{
  T *data() { return elem_.data(); }
  const T *data() const { return elem_.data(); }
  ///@}

  std::size_t n_rows() const { return desc_.extents[0]; }
  std::size_t n_cols() const { return desc_.extents[1]; }

  T &operator()(std::size_t i, std::size_t j) {
    if (!desc_.other_half(i, j))
      return *(data() + desc_(i, j));
    else
      return elem_in_other_half(i, j);
  }

  const T &operator()(std::size_t i, std::size_t j) const {
    if (!desc_.other_half(i, j))
      return *(data() + desc_(i, j));
    else
      return elem_in_other_half(i, j);
  }

 protected:
  std::vector<T> elem_;
  TRI desc_;

  virtual T &elem_in_other_half(std::size_t i, std::size_t j) = 0;
  virtual const T &elem_in_other_half(std::size_t i, std::size_t j) const = 0;
};

template<typename T, typename TRI>
class SymmetricMatrix : public PackedMatrix<T, TRI> {
 public:
  SymmetricMatrix(std::size_t n) : PackedMatrix<T, TRI>{n} {}

  T &elem_in_other_half(std::size_t i, std::size_t j) {
    return *(this->data() + this->desc_(j, i));
  }

  const T &elem_in_other_half(std::size_t i, std::size_t j) const {
    return *(this->data() + this->desc_(j, i));
  }
};

template<typename T, typename TRI>
class TriangularMatrix : public PackedMatrix<T, TRI> {
 public:
  TriangularMatrix(std::size_t n) : PackedMatrix<T, TRI>{n} {}

  T &elem_in_other_half(std::size_t i, std::size_t j) {
    return 0;
  }

  const T &elem_in_other_half(std::size_t i, std::size_t j) const {
    return 0;
  }
};

template<typename T, typename TRI>
std::ostream &operator<<(std::ostream &os, const PackedMatrix<T, TRI> &m) {
  for (std::size_t i = 0; i != m.n_rows(); ++i) {
    for (std::size_t j = 0; j != m.n_cols(); ++j) {
      os << m(i, j) << "\t";
    }
    os << std::endl;
  }

  return os << std::endl;
}

#endif // SLAB_MATRIX_PACKED_MATRIX_H
