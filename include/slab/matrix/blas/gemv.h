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
/// @brief Matrix-vector product using a general matrix

#ifndef SLAB_MATRIX_BLAS_GEMV_H_
#define SLAB_MATRIX_BLAS_GEMV_H_

namespace slab {

/// @addtogroup blas_interface BLAS INTERFACE
/// @{

/// @addtogroup blas_level2 BLAS Level 2
/// @{

template <typename T>
inline void blas_gemv(const CBLAS_TRANSPOSE trans, const T &alpha,
                      const MatrixBase<T, 2> &a, const MatrixBase<T, 1> &x,
                      const T &beta, MatrixBase<T, 1> &y) {
  const int m = y.n_rows();
  const int n = x.n_rows();

  const int lda = a.n_cols();

  const int incx = x.descriptor().strides[0];
  const int incy = y.descriptor().strides[0];

  if (is_double<T>::value) {
    cblas_dgemv(CblasRowMajor, trans, m, n, (const double)alpha,
                (const double *)(a.data() + a.descriptor().start), lda,
                (const double *)(x.data() + x.descriptor().start), incx,
                (const double)beta, (double *)(y.data() + y.descriptor().start),
                incy);
  } else if (is_float<T>::value) {
    cblas_sgemv(CblasRowMajor, trans, m, n, (const float)alpha,
                (const float *)(a.data() + a.descriptor().start), lda,
                (const float *)(x.data() + x.descriptor().start), incx,
                (const float)beta, (float *)(y.data() + y.descriptor().start),
                incy);
  } else if (is_complex_double<T>::value) {
    cblas_zgemv(
        CblasRowMajor, trans, m, n, reinterpret_cast<const double *>(&alpha),
        reinterpret_cast<const double *>(a.data() + a.descriptor().start), lda,
        reinterpret_cast<const double *>(x.data() + x.descriptor().start), incx,
        reinterpret_cast<const double *>(&beta),
        reinterpret_cast<double *>(y.data() + y.descriptor().start), incy);
  } else if (is_complex_float<T>::value) {
    cblas_cgemv(CblasRowMajor, trans, m, n, reinterpret_cast<const float *>(&alpha),
                reinterpret_cast<const float *>(a.data() + a.descriptor().start),
                lda,
                reinterpret_cast<const float *>(x.data() + x.descriptor().start),
                incx, reinterpret_cast<const float *>(&beta),
                reinterpret_cast<float *>(y.data() + y.descriptor().start), incy);
  } else {
    err_quit("blas_gemv(): unsupported element type.");
  }
}
 
/// @}
/// @} BLAS INTERFACE

}  // namespace slab
 
#endif
