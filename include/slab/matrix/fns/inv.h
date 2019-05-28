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

/// @file inv.h
/// @brief inverse of general square matrix.

#ifndef _SLAB_MATRIX_FNS_INV_H
#define _SLAB_MATRIX_FNS_INV_H

_SLAB_BEGIN_NAMESPACE

#ifndef _SLAB_USE_NO_LAPACK
template <typename T>
inline bool inverse(Matrix<T, 2> &b, const Matrix<T, 2> &a) {
  assert(a.n_rows() == a.n_cols());

  int info;
  const int m = a.n_rows();
  const int n = a.n_cols();
  const int lda = a.n_cols();
  Matrix<int, 1> ipiv(n);

  b = a;
  if (is_double<T>::value) {
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, (double *)b.data(), lda,
                          ipiv.data());
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, (double *)b.data(), lda, ipiv.data());
  } else if (is_float<T>::value) {
    info = LAPACKE_sgetrf(LAPACK_ROW_MAJOR, m, n, (float *)b.data(), lda,
                          ipiv.data());
    LAPACKE_sgetri(LAPACK_ROW_MAJOR, n, (float *)b.data(), lda, ipiv.data());
  } else if (is_complex_double<T>::value) {
    info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m, n,
                          (lapack_complex_double *)b.data(), lda, ipiv.data());
    LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, (lapack_complex_double *)b.data(), lda,
                   ipiv.data());
  } else if (is_complex_float<T>::value) {
    info = LAPACKE_cgetrf(LAPACK_ROW_MAJOR, m, n,
                          (lapack_complex_float *)b.data(), lda, ipiv.data());
    LAPACKE_cgetri(LAPACK_ROW_MAJOR, n, (lapack_complex_float *)b.data(), lda,
                   ipiv.data());
  } else {
    _SLAB_ERROR("inverse(): unspported element type.");
  }

  if (info) return false;

  return true;
}

template <typename T>
inline bool inverse(Matrix<T, 2> &b, const MatrixRef<T, 2> &a) {
  assert(a.n_rows() == a.n_cols());

  int info;
  const int m = a.n_rows();
  const int n = a.n_cols();
  const int lda = a.n_cols();
  Matrix<int, 1> ipiv(n);

  b = a;
  if (is_double<T>::value) {
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, (double *)b.data(), lda,
                          ipiv.data());
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, (double *)b.data(), lda, ipiv.data());
  } else if (is_float<T>::value) {
    info = LAPACKE_sgetrf(LAPACK_ROW_MAJOR, m, n, (float *)b.data(), lda,
                          ipiv.data());
    LAPACKE_sgetri(LAPACK_ROW_MAJOR, n, (float *)b.data(), lda, ipiv.data());
  } else if (is_complex_double<T>::value) {
    info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m, n,
                          (lapack_complex_double *)b.data(), lda, ipiv.data());
    LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, (lapack_complex_double *)b.data(), lda,
                   ipiv.data());
  } else if (is_complex_float<T>::value) {
    info = LAPACKE_cgetrf(LAPACK_ROW_MAJOR, m, n,
                          (lapack_complex_float *)b.data(), lda, ipiv.data());
    LAPACKE_cgetri(LAPACK_ROW_MAJOR, n, (lapack_complex_float *)b.data(), lda,
                   ipiv.data());
  } else {
    _SLAB_ERROR("inverse(): unspported element type.");
  }

  if (info) return false;

  return true;
}

#else
template <typename T>
void lu(Matrix<T, 2> &L, Matrix<T, 2> &U, const Matrix<T, 2> &X);

template <typename T, typename T1>
double dot_product(const MatrixRef<T, 1> &a, const MatrixRef<T1, 1> &b) {
  return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

template <typename T>
Matrix<T, 1> back_substitution(const Matrix<T, 2> &A, const Matrix<T, 1> &b) {
  const std::size_t n = A.n_rows();
  Matrix<T, 1> x(n);

  for (int i = n - 1; i >= 0; --i) {
    T s = b(i) - dot_product(A[i](slice(i + 1)), x(slice(i + 1)));
    if (T m = A(i, i))
      x(i) = s / m;
    else
      ;
  }
  return x;
}

template <typename T>
Matrix<T, 1> forward_substitution(const Matrix<T, 2> &A,
                                  const Matrix<T, 1> &b) {
  const std::size_t n = A.n_rows();
  Matrix<T, 1> x(n);

  for (int i = 0; i < n; ++i) {
    T s = {};
    if (i == 0)
      s = b(i);
    else {
      s = b(i) - dot_product(A[i](slice(0, i)), x(slice(0, i)));
    }
    if (T m = A(i, i))
      x(i) = s / m;
    else
      ;
  }
  return x;
}

template <typename T>
inline bool inverse(Matrix<T, 2> &b, const Matrix<T, 2> &a) {
  std::size_t n = a.n_rows();
  b = zeros<Matrix<T, 2> >(n, n);

  Matrix<T, 2> L, U;
  lu(L, U, a);

  for (std::size_t i = 0; i != n; ++i) {
    Matrix<T, 1> vec_one = zeros<Matrix<T, 1> >(n);
    vec_one(i) = T{1};

    Matrix<T, 1> tmp;
    tmp = forward_substitution(L, vec_one);
    b.col(i) = back_substitution(U, tmp);
  }

  return true;
}

#endif

template <typename T>
inline Matrix<T, 2> inverse(const Matrix<T, 2> &a) {
  assert(a.n_rows() == a.n_cols());

  Matrix<T, 2> res;
  inverse(res, a);

  return res;
}

template <typename T>
inline Matrix<T, 2> inverse(const MatrixRef<T, 2> &a) {
  assert(a.n_rows() == a.n_cols());

  Matrix<T, 2> res;
  inverse(res, a);

  return res;
}

template <typename T>
inline bool inv(Matrix<T, 2> &b, const Matrix<T, 2> &a) {
  return inverse(b, a);
}

template <typename T>
inline bool inv(Matrix<T, 2> &b, const MatrixRef<T, 2> &a) {
  return inverse(b, a);
}

template <typename T>
inline Matrix<T, 2> inv(const Matrix<T, 2> &a) {
  return inverse(a);
}

template <typename T>
inline Matrix<T, 2> inv(const MatrixRef<T, 2> &a) {
  return inverse(a);
}

_SLAB_END_NAMESPACE

#endif  // _SLAB_MATRIX_FNS_INV_H
