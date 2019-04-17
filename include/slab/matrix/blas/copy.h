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

/// @file copy.h
/// @brief C++ template wrapper for C functions cblas_?copy

#ifndef _SLAB_MATRIX_BLAS_COPY_H
#define _SLAB_MATRIX_BLAS_COPY_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup blas_interface BLAS Interface
/// @{

/// @addtogroup blas_level1 BLAS Level 1
/// @{

/// @brief Copies vector to another vector.
///
/// The copy routines perform a vector-vector operation defined as
/// \f[
/// x = y
/// \f]
/// where \f$x\f$ and \f$y\f$ are vectors.
///
/// @param x Vector with type vec/fvec/cx_vec/cx_fvec.
/// @param y Vector with type vec/fvec/cx_vec/cx_fvec.
/// @return Void.
///
template <typename T>
inline void blas_copy(const MatrixBase<T, 1> &x, Matrix<T, 1> &y) {
  if (std::addressof(x) == std::addressof(y)) return;
  y.clear();
  y = Matrix<T, 1>(x.size());

  const std::size_t n = x.size();
  const std::size_t incx = x.descriptor().strides[0];
  const std::size_t incy = y.descriptor().strides[0];
  const T *x_ptr = x.data() + x.descriptor().start;
  T *y_ptr = y.data() + y.descriptor().start;

  if (is_double<T>::value) {
    cblas_dcopy((const SLAB_INT)n, (const double *)x_ptr, (const SLAB_INT)incx,
                (double *)y_ptr, (const SLAB_INT)incy);
  } else if (is_complex_double<T>::value) {
    cblas_zcopy((const SLAB_INT)n, (const void *)x_ptr, (const SLAB_INT)incx,
                (void *)y_ptr, (const SLAB_INT)incy);
  }
#ifndef _SLAB_USE_R_BLAS
  else if (is_float<T>::value) {
    cblas_scopy((const SLAB_INT)n, (const float *)x_ptr, (const SLAB_INT)incx,
                (float *)y_ptr, (const SLAB_INT)incy);
  } else if (is_complex_float<T>::value) {
    cblas_ccopy((const SLAB_INT)n, (const void *)x_ptr, (const SLAB_INT)incx,
                (void *)y_ptr, (const SLAB_INT)incy);
  }
#endif
  else {
    _SLAB_ERROR("blas_copy(): unsupported element type.");
  }
}

/// @} BLAS Level 1
/// @} BLAS Interface

_SLAB_END_NAMESPACE

#endif
