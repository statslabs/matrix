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
/// @brief C++ template wrapper for C functions cblas_?gemm

#ifndef SLAB_MATRIX_BLAS_GEMM_H_
#define SLAB_MATRIX_BLAS_GEMM_H_

namespace slab {

/// @addtogroup blas_interface BLAS Interface
/// @{

/// @addtogroup blas_level3 BLAS Level 3
/// @{

/// @brief Computes a matrix-matrix product with general matrices.
///
/// The gemm routines compute a scalar-matrix-matrix product and add the result
/// to a scalar-matrix product, with general matrices. The operation is defined
/// as
/// \f[
/// C := alpha*op(A)*op(B) + beta
/// \f]
/// where \f$op(X)\f$ is one of \f$op(X) = X\f$, or \f$op(X) = X^T\f$, or
/// \f$op(X) = X^H\f$, \f$alpha\f$ and \f$beta\f$ are scalars, \f$A\f$, \f$B\f$
/// and \f$C\f$ are matrices: \f$op(A)\f$ is an m-by-k matrix, \f$op(B)\f$ is a
/// k-by-n matrix, \f$C\f$ is an m-by-n matrix.
///
template <typename T, typename T1, typename T2>
inline void blas_gemm(const CBLAS_TRANSPOSE transa,
                      const CBLAS_TRANSPOSE transb, const T1 &alpha,
                      const MatrixBase<T, 2> &a, const MatrixBase<T, 2> &b,
                      const T2 &beta, MatrixBase<T, 2> &c) {
  static_assert(Convertible<T1, T>(),
                "blas_gemm(): incompatible element type for alpha");
  static_assert(Convertible<T2, T>(),
                "blas_gemm(): incompatible element type for beta");

  const std::size_t m = c.n_rows();
  const std::size_t n = c.n_cols();
  std::size_t k = a.n_cols();

  if (transa != CblasNoTrans) k = a.n_rows();

  const std::size_t lda = a.n_cols();
  const std::size_t ldb = b.n_cols();
  const std::size_t ldc = c.n_cols();

  if (is_double<T>::value) {
    cblas_dgemm(
        CblasRowMajor, transa, transb, (const int)m, (const int)n, (const int)k,
        (const double)alpha, (const double *)(a.data() + a.descriptor().start),
        (const int)lda, (const double *)(b.data() + b.descriptor().start),
        (const int)ldb, (const double)beta,
        (double *)(c.data() + c.descriptor().start), (const int)ldc);
  } else if (is_complex_double<T>::value) {
    cblas_zgemm(
        CblasRowMajor, transa, transb, (const int)m, (const int)n, (const int)k,
        reinterpret_cast<const double *>(&alpha),
        reinterpret_cast<const double *>(a.data() + a.descriptor().start),
        (const int)lda,
        reinterpret_cast<const double *>(b.data() + b.descriptor().start),
        (const int)ldb, reinterpret_cast<const double *>(&beta),
        reinterpret_cast<double *>(c.data() + c.descriptor().start),
        (const int)ldc);
  }
#ifndef USE_R_BLAS
  else if (is_float<T>::value) {
    cblas_sgemm(
        CblasRowMajor, transa, transb, (const int)m, (const int)n, (const int)k,
        (const float)alpha, (const float *)(a.data() + a.descriptor().start),
        (const int)lda, (const float *)(b.data() + b.descriptor().start),
        (const int)ldb, (const float)beta,
        (float *)(c.data() + c.descriptor().start), (const int)ldc);

  } else if (is_complex_float<T>::value) {
    cblas_cgemm(
        CblasRowMajor, transa, transb, (const int)m, (const int)n, (const int)k,
        reinterpret_cast<const float *>(&alpha),
        reinterpret_cast<const float *>(a.data() + a.descriptor().start),
        (const int)lda,
        reinterpret_cast<const float *>(b.data() + b.descriptor().start),
        (const int)ldb, reinterpret_cast<const float *>(&beta),
        reinterpret_cast<float *>(c.data() + c.descriptor().start),
        (const int)ldc);
  }
#endif
  else {
    SLAB_ERROR("blas_gemm(): unsupported element type.");
  }
}

/// @}
/// @} BLAS Interface

}  // namespace slab

#endif
