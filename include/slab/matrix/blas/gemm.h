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

/// @file gemm.h
/// @brief Computes a matrix-matrix product with general matrices.

#ifndef SLAB_MATRIX_BLAS_GEMM_H_
#define SLAB_MATRIX_BLAS_GEMM_H_

namespace slab {

/// @addtogroup blas_interface BLAS INTERFACE
/// @{

/// @addtogroup blas_level3 BLAS Level 3
/// @{

template <typename T, typename T1, typename T2>
inline void blas_gemm(const CBLAS_TRANSPOSE transa,
                      const CBLAS_TRANSPOSE transb, const T1 &alpha,
                      const MatrixBase<T, 2> &a, const MatrixBase<T, 2> &b,
                      const T2 &beta, MatrixBase<T, 2> &c) {
  static_assert(Convertible<T1, T>(),
                "blas_gemm(): incompatible element type for alpha");
  static_assert(Convertible<T2, T>(),
                "blas_gemm(): incompatible element type for beta");

  const int m = c.n_rows();
  const int n = c.n_cols();
  int k = a.n_cols();

  if (transa != CblasNoTrans) k = a.n_rows();

  const int lda = a.n_cols();
  const int ldb = b.n_cols();
  const int ldc = c.n_cols();

  if (is_double<T>::value) {
    cblas_dgemm(CblasRowMajor, transa, transb, m, n, k, (const double)alpha,
                (const double *)(a.data() + a.descriptor().start), lda,
                (const double *)(b.data() + b.descriptor().start), ldb,
                (const double)beta, (double *)(c.data() + c.descriptor().start),
                ldc);
  } else if (is_float<T>::value) {
    cblas_sgemm(CblasRowMajor, transa, transb, m, n, k, (const float)alpha,
                (const float *)(a.data() + a.descriptor().start), lda,
                (const float *)(b.data() + b.descriptor().start), ldb,
                (const float)beta, (float *)(c.data() + c.descriptor().start),
                ldc);
  } else if (is_complex_double<T>::value) {
    cblas_zgemm(
        CblasRowMajor, transa, transb, m, n, k,
        reinterpret_cast<const double *>(&alpha),
        reinterpret_cast<const double *>(a.data() + a.descriptor().start), lda,
        reinterpret_cast<const double *>(b.data() + b.descriptor().start), ldb,
        reinterpret_cast<const double *>(&beta),
        reinterpret_cast<double *>(c.data() + c.descriptor().start), ldc);
  } else if (is_complex_float<T>::value) {
    cblas_cgemm(
        CblasRowMajor, transa, transb, m, n, k,
        reinterpret_cast<const float *>(&alpha),
        reinterpret_cast<const float *>(a.data() + a.descriptor().start), lda,
        reinterpret_cast<const float *>(b.data() + b.descriptor().start), ldb,
        reinterpret_cast<const float *>(&beta),
        reinterpret_cast<float *>(c.data() + c.descriptor().start), ldc);
  } else {
    err_quit("blas_gemm(): unsupported element type.");
  }
}

/// @}
/// @} BLAS INTERFACE

}  // namespace slab

#endif
