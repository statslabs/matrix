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

/// @file rot.h
/// @brief C++ template wrapper for C functions cblas_?rot

#ifndef _SLAB_MATRIX_BLAS_ROT_H
#define _SLAB_MATRIX_BLAS_ROT_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup blas_interface BLAS Interface
/// @{

/// @addtogroup blas_level1 BLAS Level 1
/// @{

/// @brief Performs rotation of points in the plane.
///
/// Given two complex vectors x and y, each vector element of these vectors is
/// replaced as follows: \f[ x_i = c*x_i + s*y_i \f] and \f[ y_i = c*y_i - s*x_i
/// \f]
/// where \f$x\f$ is a vector with n elements.
///
/// @param x Vector with type vec/fvec.
/// @param y Vector with type vec/fvec.
/// @param c A scalar.
/// @param s A scalar.
/// @return Void.
///
template <typename T>
inline void blas_rot(Matrix<T, 1> &x, Matrix<T, 1> &y, const T c, const T s) {
  _SLAB_ASSERT(x.size() == y.size(),
               "blas_rot(): incompatible vector dimensions");

  const std::size_t n = x.size();
  const std::size_t incx = x.descriptor().strides[0];
  const std::size_t incy = y.descriptor().strides[0];
  const T *x_ptr = x.data() + x.descriptor().start;
  const T *y_ptr = y.data() + y.descriptor().start;

  if (is_double<T>::value) {
    cblas_drot((const SLAB_INT)n, (double *)x_ptr, (const SLAB_INT)incx,
               (double *)y_ptr, (const SLAB_INT)incy, (const double)c,
               (const double)s);
  }
#ifndef _SLAB_USE_R_BLAS
  else if (is_float<T>::value) {
    cblas_srot((const SLAB_INT)n, (float *)x_ptr, (const SLAB_INT)incx,
               (float *)y_ptr, (const SLAB_INT)incy, (const float)c,
               (const float)s);
  }
#endif
  else {
    _SLAB_ERROR("blas_rot(): unsupported element type.");
  }
}

#ifndef _SLAB_USE_R_BLAS

/// @brief Performs rotation of points in the plane.
///
/// Given two complex vectors x and y, each vector element of these vectors is
/// replaced as follows: \f[ x_i = c*x_i + s*y_i \f] and \f[ y_i = c*y_i - s*x_i
/// \f]
/// where \f$x\f$ is a vector with n elements.
///
/// @param x Vector with type cx_vec/cx_fvec.
/// @param y Vector with type cx_vec/cx_fvec.
/// @param c A scalar.
/// @param s A scalar.
/// @return Void.
///
template <typename T>
inline void blas_rot(Matrix<std::complex<T>, 1> &x,
                     Matrix<std::complex<T>, 1> &y, const T c, const T s) {
  _SLAB_ASSERT(x.size() == y.size(),
               "blas_rot(): incompatible vector dimensions");

  const std::size_t n = x.size();
  const std::size_t incx = x.descriptor().strides[0];
  const std::size_t incy = y.descriptor().strides[0];
  const std::complex<T> *x_ptr = x.data() + x.descriptor().start;
  const std::complex<T> *y_ptr = y.data() + y.descriptor().start;

  if (is_double<T>::value) {
    cblas_zdrot((const SLAB_INT)n, (void *)x_ptr, (const SLAB_INT)incx,
                (void *)y_ptr, (const SLAB_INT)incy, (const double)c,
                (const double)s);
  } else if (is_float<T>::value) {
    cblas_csrot((const SLAB_INT)n, (void *)x_ptr, (const SLAB_INT)incx,
                (void *)y_ptr, (const SLAB_INT)incy, (const float)c,
                (const float)s);
  } else {
    _SLAB_ERROR("blas_rot(): unsupported element type.");
  }
}

#endif

/// @} BLAS Level 1
/// @} BLAS Interface

_SLAB_END_NAMESPACE

#endif
