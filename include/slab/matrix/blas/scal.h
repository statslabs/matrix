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

/// @file scal.h
/// @brief C++ template wrapper for C functions cblas_?scal

#ifndef _SLAB_MATRIX_BLAS_SCAL_H
#define _SLAB_MATRIX_BLAS_SCAL_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup blas_interface BLAS Interface
/// @{

/// @addtogroup blas_level1 BLAS Level 1
/// @{

/// @brief Computes the product of a vector by a scalar.
///
template <typename T>
inline void blas_scal(const T a, Matrix<T, 1> &x) {
  const int n = x.size();
  const int incx = x.descriptor().strides[0];

  if (is_double<T>::value) {
    cblas_dscal(n, (const double)a, (double *)(x.data() + x.descriptor().start),
                incx);
  } else if (is_complex_double<T>::value) {
    cblas_zdscal(n, (const double)a,
                 reinterpret_cast<double *>(x.data() + x.descriptor().start),
                 incx);
  }
#ifndef _SLAB_USE_R_BLAS
  else if (is_float<T>::value) {
    cblas_sscal(n, (const float)a, (float *)(x.data() + x.descriptor().start),
                incx);

  } else if (is_complex_float<T>::value) {
    cblas_csscal(n, (const float)a,
                 reinterpret_cast<float *>(x.data() + x.descriptor().start),
                 incx);
  }
#endif
  else {
    _SLAB_ERROR("blas_scal(): unsupported element type.");
  }
}

/// @brief Computes the product of a vector by a scalar.
///
template <typename T>
inline void blas_scal(const std::complex<T> &a, Matrix<std::complex<T>, 1> &x) {
  const int n = x.size();
  const int incx = x.descriptor().strides[0];

  if (is_double<T>::value) {
    cblas_zscal(n, reinterpret_cast<const double *>(&a),
                reinterpret_cast<double *>(x.data() + x.descriptor().start),
                incx);
  }
#ifndef _SLAB_USE_R_BLAS
  else if (is_float<T>::value) {
    cblas_cscal(n, reinterpret_cast<const float *>(&a),
                reinterpret_cast<float *>(x.data() + x.descriptor().start),
                incx);
  }
#endif
  else {
    _SLAB_ERROR("blas_scal(): unsupported element type.");
  }
}

/// @} BLAS Level 1
/// @} BLAS Interface

_SLAB_END_NAMESPACE

#endif
