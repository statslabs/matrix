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

/// @file spr.h
/// @brief C++ template wrapper for C functions cblas_?spr

#ifndef _SLAB_MATRIX_BLAS_SPR_H
#define _SLAB_MATRIX_BLAS_SPR_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup blas_interface BLAS Interface
/// @{

/// @addtogroup blas_level2 BLAS Level 2
/// @{

/// @brief Performs a rank-1 update of a symmetric packed matrix.
///
template <typename T, typename TRI>
inline void blas_spr(const T &alpha, const MatrixBase<T, 1> &x,
                     SymmetricMatrix<T, TRI> &ap) {
  assert(x.size() == ap.n_rows());

  CBLAS_UPLO uplo;
  if (is_upper<TRI>::value)
    uplo = CblasUpper;
  else if (is_lower<TRI>::value)
    uplo = CblasLower;

  const int incx = x.descriptor().strides[0];

  if (is_double<T>::value) {
    cblas_dspr(CblasRowMajor, uplo, x.size(), (const double)alpha,
               (const double *)x.data(), incx, (double *)ap.data());
  }
#ifndef _SLAB_USE_R_BLAS
  else if (is_float<T>::value) {
    cblas_sspr(CblasRowMajor, uplo, x.size(), (const float)alpha,
               (const float *)x.data(), incx, (float *)ap.data());
  }
#endif
  else {
    _SLAB_ERROR("blas_spr(): unsupported element type.");
  }
}

/// @} BLAS Level 2
/// @} BLAS Interface

_SLAB_END_NAMESPACE

#endif
