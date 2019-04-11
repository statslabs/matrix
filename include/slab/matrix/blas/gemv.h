//
// Copyright 2019 The Statslabs Authors.
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

#ifndef SLAB_MATRIX_BLAS_GEMV_H_
#define SLAB_MATRIX_BLAS_GEMV_H_

namespace slab {

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
template <typename T>
inline void blas_gemv(const CBLAS_TRANSPOSE trans, const T &alpha,
                      const MatrixBase<T, 2> &a, const MatrixBase<T, 1> &x,
                      const T &beta, MatrixBase<T, 1> &y) {
  const std::size_t m = y.n_rows();
  const std::size_t n = x.n_rows();

  const std::size_t lda = a.n_cols();

  const std::size_t incx = x.descriptor().strides[0];
  const std::size_t incy = y.descriptor().strides[0];

  if (is_double<T>::value) {
    cblas_dgemv(
        CblasRowMajor, trans, (const int)m, (const int)n, (const double)alpha,
        (const double *)(a.data() + a.descriptor().start), (const int)lda,
        (const double *)(x.data() + x.descriptor().start), (const int)incx,
        (const double)beta, (double *)(y.data() + y.descriptor().start),
        (const int)incy);
  } else if (is_complex_double<T>::value) {
    cblas_zgemv(
        CblasRowMajor, trans, (const int)m, (const int)n,
        reinterpret_cast<const double *>(&alpha),
        reinterpret_cast<const double *>(a.data() + a.descriptor().start),
        (const int)lda,
        reinterpret_cast<const double *>(x.data() + x.descriptor().start),
        (const int)incx, reinterpret_cast<const double *>(&beta),
        reinterpret_cast<double *>(y.data() + y.descriptor().start),
        (const int)incy);
  }
#ifndef _SLAB_USE_R_BLAS
  else if (is_float<T>::value) {
    cblas_sgemv(
        CblasRowMajor, trans, (const int)m, (const int)n, (const float)alpha,
        (const float *)(a.data() + a.descriptor().start), (const int)lda,
        (const float *)(x.data() + x.descriptor().start), (const int)incx,
        (const float)beta, (float *)(y.data() + y.descriptor().start),
        (const int)incy);
  } else if (is_complex_float<T>::value) {
    cblas_cgemv(
        CblasRowMajor, trans, (const int)m, (const int)n,
        reinterpret_cast<const float *>(&alpha),
        reinterpret_cast<const float *>(a.data() + a.descriptor().start),
        (const int)lda,
        reinterpret_cast<const float *>(x.data() + x.descriptor().start),
        (const int)incx, reinterpret_cast<const float *>(&beta),
        reinterpret_cast<float *>(y.data() + y.descriptor().start),
        (const int)incy);
  }
#endif
  else {
    _SLAB_ERROR("blas_gemv(): unsupported element type.");
  }
}

/// @}
/// @} BLAS Interface

}  // namespace slab

#endif
