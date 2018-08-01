#ifndef SLAB_MATRIX_PACKED_MATRIX_H_
#define SLAB_MATRIX_PACKED_MATRIX_H_

// Triangular matrix type
struct lower_tag {};
struct upper_tag {};
struct unit_lower_tag : public lower_tag {};
struct unit_upper_tag : public upper_tag {};

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

  PackedMatrix(std::size_t);

  //! "flat" element access
  ///@{
  T *data() { return elem_.data(); }
  const T *data() { return elem_.data(); }
  ///@}

  T &operator() (std::size_t i, std::size_t j) {
    if (!desc_.other_half(i, j))
      return *(data() + desc_(i, j, n_rows()));
    else
      return *(data() + desc_(j, i, n_rows()));
  }

  const T &operator() (std::size_t i, std::size_t j) const {
    if (!desc_.other_half(i, j))
      return *(data() + desc_(i, j, n_rows()));
    else
      return *(data() + desc_(j, i, n_rows()));
  }

 protected:
  std::vector<T> elem_;
  TRI desc_;
};

#endif // SLAB_MATRIX_PACKED_MATRIX_H
