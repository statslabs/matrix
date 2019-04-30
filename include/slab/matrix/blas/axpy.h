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

/// @file axpy.h
/// @brief C++ template wrapper for C functions cblas_?axpy

#ifndef _SLAB_MATRIX_BLAS_AXPY_H
#define _SLAB_MATRIX_BLAS_AXPY_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup blas_interface BLAS Interface
/// @{

/// @addtogroup blas_level1 BLAS Level 1
/// @{

/// @brief Computes a vector-scalar product and adds the result to a vector.
///
/// The axpy routines perform a vector-vector operation defined as
/// \f[
/// y := a*x + y
/// \f]
/// where \f$a\f$ is a scalar, \f$x\f$ and \f$y\f$ are vectors each
/// with a number of elements that equals n.
///
/// @param a Specifies the scalar a.
/// @param x Vector with type vec/fvec/cx_vec/cx_fvec.
/// @param y Vector with typeaxp vec/fvec/cx_vec/cx_fvec.
/// @return Void.
///
template <typename T, typename T1>
inline void blas_axpy(const T1 &a, const MatrixBase<T, 1> &x,
                      MatrixBase<T, 1> &y) {
  static_assert(Convertible<T1, T>(),
                "blas_axpy(): incompatible element type for alpha");
  _SLAB_ASSERT(x.size() == y.size(),
               "blas_axpy(): incompatible vector dimensions");
  _SLAB_ASSERT(std::addressof(x) != std::addressof(y),
               "blas_axpy(): x and y should not be the same object");

  const std::size_t n = x.size();
  const std::size_t incx = x.descriptor().strides[0];
  const std::size_t incy = y.descriptor().strides[0];
  const T *x_ptr = x.data() + x.descriptor().start;
  T *y_ptr = y.data() + y.descriptor().start;

  if (is_double<T>::value) {
    cblas_daxpy((const SLAB_INT)n, (const double)a, (const double *)x_ptr,
                (const SLAB_INT)incx, (double *)y_ptr, (const SLAB_INT)incy);
  } else if (is_complex_double<T>::value) {
    cblas_zaxpy((const SLAB_INT)n, (const void *)(&a), (const void *)x_ptr,
                (const SLAB_INT)incx, (void *)y_ptr, (const SLAB_INT)incy);
  }
#ifndef _SLAB_USE_R_BLAS
  else if (is_float<T>::value) {
    cblas_saxpy((const SLAB_INT)n, (const float)a, (const float *)x_ptr,
                (const SLAB_INT)incx, (float *)y_ptr, (const SLAB_INT)incy);
  } else if (is_complex_float<T>::value) {
    cblas_caxpy((const SLAB_INT)n, (const void *)(&a), (const void *)x_ptr,
                (const SLAB_INT)incx, (void *)y_ptr, (const SLAB_INT)incy);
  }
#endif
  else {
    _SLAB_ERROR("blas_axpy(): unsupported element type.");
  }
}

/// @} BLAS Level 1
/// @} BLAS Interface

_SLAB_END_NAMESPACE

#endif
