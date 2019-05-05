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

/// @file rotm.h
/// @brief C++ template wrapper for C functions cblas_?rotm

#ifndef _SLAB_MATRIX_BLAS_ROTM_H
#define _SLAB_MATRIX_BLAS_ROTM_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup blas_interface BLAS Interface
/// @{

/// @addtogroup blas_level1 BLAS Level 1
/// @{

/// @brief Performs modified Givens rotation of points in the plane.
///
/// Given two vectors x and y, each vector element of these vectors is replaced
/// as follows:
///
/// for i=1 to n, where H is a modified Givens transformation matrix whose
/// values are stored in the param[1] through param[4] array. See discussion on
/// the param argument.
///
/// @param x Vector with type vec/fvec.
/// @param y Vector with type vec/fvec.
/// @param param Vector with type vec/fvec, size 5.
/// @return Void.
///
template <typename T>
inline void blas_rotm(Matrix<T, 1> &x, Matrix<T, 1> &y, Matrix<T, 1> &param) {
  _SLAB_ASSERT(x.size() == y.size(),
               "blas_rotm(): incompatible vector dimensions");
  _SLAB_ASSERT(param.size() == 5,
               "blas_rotm(): param should be a vector of size 5");

  const std::size_t n = x.size();
  const std::size_t incx = x.descriptor().strides[0];
  const std::size_t incy = y.descriptor().strides[0];
  T *x_ptr = x.data() + x.descriptor().start;
  T *y_ptr = y.data() + y.descriptor().start;
  const T *param_ptr = param.data() + param.descriptor().start;

  if (is_double<T>::value) {
    cblas_drotm((const SLAB_INT)n, (double *)x_ptr, (const SLAB_INT)incx,
                (double *)y_ptr, (const SLAB_INT)incy,
                (const double *)param_ptr);
  }
#ifndef _SLAB_USE_R_BLAS
  else if (is_float<T>::value) {
    cblas_srotm((const SLAB_INT)n, (float *)x_ptr, (const SLAB_INT)incx,
                (float *)y_ptr, (const SLAB_INT)incy, (const float *)param_ptr);
  }
#endif
  else {
    _SLAB_ERROR("blas_rotm(): unsupported element type.");
  }
}

/// @} BLAS Level 1
/// @} BLAS Interface

_SLAB_END_NAMESPACE

#endif
