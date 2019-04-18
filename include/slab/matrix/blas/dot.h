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

/// @file dot.h
/// @brief C++ template wrapper for C functions cblas_?dot

#ifndef _SLAB_MATRIX_BLAS_DOT_H
#define _SLAB_MATRIX_BLAS_DOT_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup blas_interface BLAS Interface
/// @{

/// @addtogroup blas_level1 BLAS Level 1
/// @{

/// @brief Computes a vector-vector dot product.
///
/// The dot routines perform a vector-vector reduction operation
/// defined as
/// \f[
/// res = \sum_{i=1}^{n} x_i \times y_i
/// \f]
/// where \f$x_i\f$ and \f$y_i\f$ are elements of vectors \f$x\f$ and
/// \f$y\f$.
///
/// @param x Vector with type vec/fvec/cx_vec/cx_fvec.
/// @param y Vector with type vec/fvec/cx_vec/cx_fvec.
/// @return The result of the dot product of \f$x\f$ and \f$y\f$.
///
template <typename T>
inline T blas_dot(const MatrixBase<T, 1> &x, const MatrixBase<T, 1> &y) {
  _SLAB_ASSERT(x.size() == y.size(),
               "blas_doc(): incompatible vector dimensions");

  const std::size_t n = x.size();
  const std::size_t incx = x.descriptor().strides[0];
  const std::size_t incy = y.descriptor().strides[0];
  const T *x_ptr = x.data() + x.descriptor().start;
  const T *y_ptr = y.data() + y.descriptor().start;

  T res = {};  // C++11: zero initialization
  if (is_double<T>::value) {
    res = cblas_ddot((const SLAB_INT)n, (const double *)x_ptr,
                     (const SLAB_INT)incx, (const double *)y_ptr,
                     (const SLAB_INT)incy);
  }
#ifndef _SLAB_USE_R_BLAS
  else if (is_float<T>::value) {
    res = cblas_sdot((const SLAB_INT)n, (const float *)x_ptr,
                     (const SLAB_INT)incx, (const float *)y_ptr,
                     (const SLAB_INT)incy);
  }
#endif
  else {
    _SLAB_ERROR("blas_dot(): unsupported element type.");
  }

  return res;
}

/// @} BLAS Level 1
/// @} BLAS Interface

_SLAB_END_NAMESPACE

#endif
