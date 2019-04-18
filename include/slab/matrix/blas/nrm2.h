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

/// @file nrm2.h
/// @brief C++ template wrapper for C functions cblas_?nrm2

#ifndef _SLAB_MATRIX_BLAS_NRM2_H
#define _SLAB_MATRIX_BLAS_NRM2_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup blas_interface BLAS Interface
/// @{

/// @addtogroup blas_level1 BLAS Level 1
/// @{

/// @brief Computes the Euclidean norm of a vector.
///
/// The nrm2 routines perform a vector reduction operation defined as
/// \f[
/// res = ||x||
/// \f]
/// where: \f$x\f$ is a vector, \f$res\f$ is a value containing the
/// Euclidean norm of the elements of \f$x\f$.
///
/// @param x Vector with type vec/fvec/cx_vec/cx_fvec.
/// @return The Euclidean norm of the vector x.
///
template <typename T>
inline double blas_nrm2(const Matrix<T, 1> &x) {
  const std::size_t n = x.size();
  const std::size_t incx = x.descriptor().strides[0];
  const T *x_ptr = x.data() + x.descriptor().start;

  double res = 0.0;
  if (is_double<T>::value) {
    res = cblas_dnrm2((const SLAB_INT)n, (const double *)x_ptr,
                      (const SLAB_INT)incx);
  } else if (is_complex_double<T>::value) {
    res = cblas_dznrm2((const SLAB_INT)n, (const void *)x_ptr,
                       (const SLAB_INT)incx);
  }
#ifndef _SLAB_USE_R_BLAS
  else if (is_float<T>::value) {
    return cblas_snrm2((const SLAB_INT)n, (const float *)x_ptr,
                       (const SLAB_INT)incx);
  } else if (is_complex_float<T>::value) {
    return cblas_scnrm2((const SLAB_INT)n, (const void *)x_ptr,
                        (const SLAB_INT)incx);
  }
#endif
  else {
    _SLAB_ERROR("blas_nrm2(): unsupported element type.");
  }

  return res;
}

/// @} BLAS Level 1
/// @} BLAS Interface

_SLAB_END_NAMESPACE

#endif
