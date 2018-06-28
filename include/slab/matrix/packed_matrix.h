#ifndef SLAB_MATRIX_PACKED_MATRIX_H_
#define SLAB_MATRIX_PACKED_MATRIX_H_

template<typename T>
class TriangularMatrix {
 public:
  using value_type = T;
  using iterator = typename std::vector<T>::iterator;
  using const_iterator  = typename std::vector<T>::const_iterator;

  TriangularMatrix() = default;
  TriangularMatrix(TriangularMatrix &&) = default;
  TriangularMatrix &operator=(TriangularMatrix &&) = default;
  TriangularMatrix(TriangularMatrix const &) = default;
  TriangularMatrix &operator=(TriangularMatrix const &) = default;
  ~TriangularMatrix() = default;

  TriangularMatrix(std::size_t, const std::initializer_list<T> &);
 private:
  std::vector<T> elem_;
  bool is_unit_tri_ = false;
};

template<typename T>
TriangularMatrix::TriangularMatrix(std::size_t n,
                                   const std::initializer_list<T> &init)
{

}

#endif // SLAB_MATRIX_PACKED_MATRIX_H
