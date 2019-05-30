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

/// @file asum.h
/// @brief C++ template wrapper for C functions cblas_?asum

#ifndef _SLAB_MATRIX_BLAS_ASUM_H
#define _SLAB_MATRIX_BLAS_ASUM_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup blas_interface BLAS Interface
/// @{

/// @addtogroup blas_level1 BLAS Level 1
/// @{

/// @brief Computes the sum of magnitudes of the vector elements.
///
/// The asum routine computes the sum of the magnitudes of elements of
/// a real vector, or the sum of magnitudes of the real and imaginary
/// parts of elements of a complex vector:
/// \f[
/// res = |Re~x_1| + |Im~x_1| + |Re~x_2| + |Im~x_2|+ ... + |Re~x_n| + |Im~x_n|,
/// \f]
/// where \f$x\f$ is a vector with n elements.
///
/// @param x Vector with type vec/fvec.
/// @return Contains the sum of magnitudes of real and imaginary parts
///         of all elements of the vector.
///
template <typename T>
inline T blas_asum(const MatrixBase<T, 1> &x) {
  const std::size_t n = x.size();
  const std::size_t incx = x.descriptor().strides[0];
  const T *x_ptr = x.data() + x.descriptor().start;

  T res = {};  // C++11: zero initialization
  if (is_double<T>::value) {
    res = cblas_dasum((const SLAB_INT)n, (const double *)x_ptr,
                      (const SLAB_INT)incx);
  }
#ifndef _SLAB_USE_R_BLAS
  else if (is_float<T>::value) {
    res = cblas_sasum((const SLAB_INT)n, (const float *)x_ptr,
                      (const SLAB_INT)incx);
  }
#endif
  else {
    _SLAB_ERROR("blas_asum(): unsupported element type.");
  }

  return res;
}

/// @brief Computes the sum of magnitudes of the vector elements.
///
/// The asum routine computes the sum of the magnitudes of elements of
/// a real vector, or the sum of magnitudes of the real and imaginary
/// parts of elements of a complex vector:
/// \f[
/// res = |Re~x_1| + |Im~x_1| + |Re~x_2| + |Im~x_2|+ ... + |Re~x_n| + |Im~x_n|,
/// \f]
/// where \f$x\f$ is a vector with n elements.
///
/// @param x Vector with type cx_vec/cx_fvec.
/// @return Contains the sum of magnitudes of real and imaginary parts
///         of all elements of the vector.
///
template <typename T>
inline T blas_asum(const MatrixBase<std::complex<T>, 1> &x) {
  const std::size_t n = x.size();
  const std::size_t incx = x.descriptor().strides[0];
  const std::complex<T> *x_ptr = x.data() + x.descriptor().start;

  T res = {};  // C++11: zero initialization
  if (is_double<T>::value) {
    res = cblas_dzasum((const SLAB_INT)n, (const void *)x_ptr,
                       (const SLAB_INT)incx);
  }
#ifndef _SLAB_USE_R_BLAS
  else if (is_float<T>::value) {
    res = cblas_scasum((const SLAB_INT)n, (const void *)x_ptr,
                       (const SLAB_INT)incx);
  }
#endif
  else {
    _SLAB_ERROR("blas_asum(): unsupported element type.");
  }

  return res;
}

/// @} BLAS Level 1
/// @} BLAS Interface

_SLAB_END_NAMESPACE

#endif
