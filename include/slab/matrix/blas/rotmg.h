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

/// @file rotmg.h
/// @brief C++ template wrapper for C functions cblas_?rotmg

#ifndef _SLAB_MATRIX_BLAS_ROTMG_H
#define _SLAB_MATRIX_BLAS_ROTMG_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup blas_interface BLAS Interface
/// @{

/// @addtogroup blas_level1 BLAS Level 1
/// @{

/// @brief Computes the parameters for a modified Givens rotation.
///
/// Given Cartesian coordinates (x1, y1) of an input vector, these routines
/// compute the components of a modified Givens transformation matrix H that
/// zeros the y-component of the resulting vector:
///
/// @param d1 Provides the scaling factor for the x-coordinate of the input
/// vector; Provides the first diagonal element of the updated matrix.
/// @param d2 Provides the scaling factor for the y-coordinate of the input
/// vector; Provides the second diagonal element of the updated matrix.
/// @param x1 Provides the x-coordinate of the input vector; Provides the
/// x-coordinate of the rotated vector before scaling.
/// @param y1 Provides the y-coordinate of the input vector.
/// @param param Vector with type vec/fvec, size 5.
/// @return Void.
///
template <typename T>
inline void blas_rotmg(Matrix<T, 1> &d1, Matrix<T, 1> &d2, Matrix<T, 1> &x1,
                       const T y1, Matrix<T, 1> &param) {
  _SLAB_ASSERT(param.size() == 5,
               "blas_rotmg(): param should be a vector of size 5");

  T *d1_ptr = d1.data() + d1.descriptor().start;
  T *d2_ptr = d2.data() + d2.descriptor().start;
  T *x1_ptr = x1.data() + x1.descriptor().start;
  T *param_ptr = param.data() + param.descriptor().start;

  if (is_double<T>::value) {
    cblas_drotmg((double *)d1_ptr, (double *)d2_ptr, (double *)x1_ptr,
                 (const double)y1, (double *)param_ptr);
  }
#ifndef _SLAB_USE_R_BLAS
  else if (is_float<T>::value) {
    cblas_srotmg((float *)d1_ptr, (float *)d2_ptr, (float *)x1_ptr,
                 (const float)y1, (float *)param_ptr);
  }
#endif
  else {
    _SLAB_ERROR("blas_rotmg(): unsupported element type.");
  }
}

/// @} BLAS Level 1
/// @} BLAS Interface

_SLAB_END_NAMESPACE

#endif
