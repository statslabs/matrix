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

/// @file swap.h
/// @brief C++ template wrapper for C functions cblas_?swap

#ifndef _SLAB_MATRIX_BLAS_SWAP_H
#define _SLAB_MATRIX_BLAS_SWAP_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup blas_interface BLAS Interface
/// @{

/// @addtogroup blas_level1 BLAS Level 1
/// @{

/// @brief Swaps a vector with another vector.
///
/// @param x a vector.
/// @param y another vector.
///
template <typename T>
inline void blas_swap(Matrix<T, 1> &x, Matrix<T, 1> &y) {
  assert(x.size() == y.size());

  const int n = x.size();
  const int incx = x.descriptor().strides[0];
  const int incy = y.descriptor().strides[0];

  if (is_double<T>::value) {
    cblas_dswap(n, (double *)(x.data() + x.descriptor().start), incx,
                (double *)(y.data() + y.descriptor().start), incy);
  } else if (is_complex_double<T>::value) {
    cblas_zswap(
        n, reinterpret_cast<double *>(x.data() + x.descriptor().start), incx,
        reinterpret_cast<double *>(y.data() + y.descriptor().start), incy);
  }
#ifndef _SLAB_USE_R_BLAS
  else if (is_float<T>::value) {
    cblas_sswap(n, (float *)(x.data() + x.descriptor().start), incx,
                (float *)(y.data() + y.descriptor().start), incy);
  } else if (is_complex_float<T>::value) {
    cblas_cswap(
        n, reinterpret_cast<float *>(x.data() + x.descriptor().start), incx,
        reinterpret_cast<float *>(y.data() + y.descriptor().start), incy);
  }
#endif
  else {
    _SLAB_ERROR("blas_swap(): unsupported element type.");
  }
}

/// @} BLAS Level 1
/// @} BLAS Interface

_SLAB_END_NAMESPACE

#endif
