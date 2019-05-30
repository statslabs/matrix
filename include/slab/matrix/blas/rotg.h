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

/// @file rotg.h
/// @brief C++ template wrapper for C functions cblas_?rotg

#ifndef _SLAB_MATRIX_BLAS_ROTG_H
#define _SLAB_MATRIX_BLAS_ROTG_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup blas_interface BLAS Interface
/// @{

/// @addtogroup blas_level1 BLAS Level 1
/// @{

/// @brief Computes the parameters for a Givens rotation.
///
/// Given the Cartesian coordinates (a, b) of a point, these routines return the
/// parameters c, s, r, and z associated with the Givens rotation. The
/// parameters c and s define a unitary matrix such that:
///
/// The parameter z is defined such that if |a| > |b|, z is s; otherwise if c is
/// not 0 z is 1/c; otherwise z is 1.
///
/// @return Void.
///
template <typename T>
void blas_rotg(Matrix<T, 1> &a, Matrix<T, 1> &b, Matrix<T, 1> &c,
               Matrix<T, 1> &s) {
  if (is_double<T>::value) {
    cblas_drotg((double *)a.data(), (double *)b.data(), (double *)c.data(),
                (double *)s.data());
  }
#ifndef _SLAB_USE_R_BLAS
  else if (is_float<T>::value) {
    cblas_srotg((float *)a.data(), (float *)b.data(), (float *)c.data(),
                (float *)s.data());
  }
#endif
  else {
    _SLAB_ERROR("blas_rotg(): unsupported element type.");
  }
}

#ifndef _SLAB_USE_R_BLAS

template <typename T>
void blas_rotg(Matrix<std::complex<T>, 1> &a,
               const Matrix<std::complex<T>, 1> &b, Matrix<T, 1> &&c,
               Matrix<std::complex<T>, 1> &s) {
  if (is_double<T>::value) {
    cblas_zrotg((void *)a.data(), (const void *)b.data(), (double *)c.data(),
                (void *)s.data());
  } else if (is_float<T>::value) {
    cblas_crotg((void *)a.data(), (const void *)b.data(), (float *)c.data(),
                (void *)s.data());
  } else {
    _SLAB_ERROR("blas_rotg(): unsupported element type.");
  }
}

#endif

/// @} BLAS Level 1
/// @} BLAS Interface

_SLAB_END_NAMESPACE

#endif
