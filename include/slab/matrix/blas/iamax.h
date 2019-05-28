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
/// @brief C++ template wrapper for C functions cblas_i?amax

#ifndef _SLAB_MATRIX_BLAS_IAMAX_H
#define _SLAB_MATRIX_BLAS_IAMAX_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup blas_interface BLAS Interface
/// @{

/// @addtogroup blas_level1 BLAS Level 1
/// @{

/// @brief Finds the index of the element with maximum absolute value.
///
template <typename T>
inline std::size_t blas_iamax(const Matrix<T, 1> &x) {
  std::size_t res = 0;
  std::size_t incx = x.descriptor().strides[0];

  if (is_double<T>::value) {
    res = cblas_idamax(x.size(),
                       (const double *)(x.data() + x.descriptor().start), incx);
  } else if (is_complex_double<T>::value) {
    res = cblas_izamax(
        x.size(),
        reinterpret_cast<const double *>(x.data() + x.descriptor().start),
        incx);
  }
#ifndef _SLAB_R_USE_BLAS
  else if (is_float<T>::value) {
    res = cblas_isamax(x.size(),
                       (const float *)(x.data() + x.descriptor().start), incx);
  } else if (is_complex_float<T>::value) {
    res = cblas_icamax(
        x.size(),
        reinterpret_cast<const float *>(x.data() + x.descriptor().start), incx);
  }
#endif
  else {
    _SLAB_ERROR("blas_iamax(): unsupported element type.");
  }

  return res;
}

/// @} BLAS Level 1
/// @} BLAS Interface

_SLAB_END_NAMESPACE

#endif
