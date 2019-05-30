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

/// @file gemv.h
/// @brief C++ template wrapper for C functions cblas_?gemv

#ifndef _SLAB_MATRIX_BLAS_GEMV_H
#define _SLAB_MATRIX_BLAS_GEMV_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup blas_interface BLAS Interface
/// @{

/// @addtogroup blas_level2 BLAS Level 2
/// @{

/// @brief Matrix-vector product using a general matrix.
///
/// The gemv routines perform a matrix-vector operation defined as:
/// \f[
/// y := alpha*A*x + beta*y
/// \f],
/// or
/// \f[
/// y := alpha*A'*x + beta*y
/// \f],
/// or
/// \f[
/// y := alpha*conjg(A')*x + beta*y
/// \f],
/// where \f$alpha\f$ and \f$beta\f$ are scalars, \f$x\f$ and \f$y\f$ are
/// vectors, \f$A\f$ is an m-by-n matrix.
///
template <typename T, typename T1, typename T2>
inline void blas_gemv(const CBLAS_TRANSPOSE trans, const T1 &alpha,
                      const MatrixBase<T, 2> &a, const MatrixBase<T, 1> &x,
                      const T2 &beta, MatrixBase<T, 1> &y) {
  static_assert(Convertible<T1, T>(),
                "blas_gemv(): incompatible element type for alpha");
  static_assert(Convertible<T2, T>(),
                "blas_gemv(): incompatible element type for beta");

  const std::size_t m = y.n_rows();
  const std::size_t n = x.n_rows();

  const std::size_t lda = a.n_cols();

  const std::size_t incx = x.descriptor().strides[0];
  const std::size_t incy = y.descriptor().strides[0];

  const T *a_ptr = a.data() + a.descriptor().start;
  const T *x_ptr = x.data() + x.descriptor().start;
  T *y_ptr = y.data() + y.descriptor().start;

  if (is_double<T>::value) {
    cblas_dgemv(CblasRowMajor, trans, (const SLAB_INT)m, (const SLAB_INT)n,
                (const double)alpha, (const double *)a_ptr, (const SLAB_INT)lda,
                (const double *)x_ptr, (const SLAB_INT)incx, (const double)beta,
                (double *)y_ptr, (const SLAB_INT)incy);
  } else if (is_complex_double<T>::value) {
    cblas_zgemv(CblasRowMajor, trans, (const SLAB_INT)m, (const SLAB_INT)n,
                (const void *)(&alpha), (const void *)a_ptr,
                (const SLAB_INT)lda, (const void *)x_ptr, (const SLAB_INT)incx,
                (const void *)(&beta), (void *)y_ptr, (const SLAB_INT)incy);
  }
#ifndef _SLAB_USE_R_BLAS
  else if (is_float<T>::value) {
    cblas_sgemv(CblasRowMajor, trans, (const SLAB_INT)m, (const SLAB_INT)n,
                (const float)alpha, (const float *)a_ptr, (const SLAB_INT)lda,
                (const float *)x_ptr, (const SLAB_INT)incx, (const float)beta,
                (float *)y_ptr, (const SLAB_INT)incy);
  } else if (is_complex_float<T>::value) {
    cblas_cgemv(CblasRowMajor, trans, (const SLAB_INT)m, (const SLAB_INT)n,
                (const void *)(&alpha), (const void *)a_ptr,
                (const SLAB_INT)lda, (const void *)x_ptr, (const SLAB_INT)incx,
                (const void *)(&beta), (void *)y_ptr, (const SLAB_INT)incy);
  }
#endif
  else {
    _SLAB_ERROR("blas_ gemv(): unsupported element type.");
  }
}

/// @}
/// @} BLAS Interface

_SLAB_END_NAMESPACE

#endif
