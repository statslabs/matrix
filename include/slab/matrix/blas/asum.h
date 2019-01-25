//
// Copyright 2019 The Statslabs Authors.
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

#ifndef SLAB_MATRIX_BLAS_ASUM_H_
#define SLAB_MATRIX_BLAS_ASUM_H_

namespace slab {

/// @addtogroup blas_interface BLAS INTERFACE
/// @{

/// @addtogroup blas_level1 BLAS Level 1
/// @{

/// @brief Computes the sum of magnitudes of the vector elements.
///
/// The asum routine computes the sum of the magnitudes of elements of
/// a real vector, or the sum of magnitudes of the real and imaginary
/// parts of elements of a complex vector:
///
/// res = |Re x1| + |Im x1| + |Re  x2| + |Im  x2|+ ... + |Re  xn| + |Im xn|,
///
/// where x is a vector with n elements.
///
/// @param x Vector.
///
/// @return Contains the sum of magnitudes of real and imaginary parts
///         of all elements of the vector.
///
template <typename T>
inline T blas_asum(const Matrix<T, 1> &x) {
  const int n = x.size();
  const int incx = x.descriptor().strides[0];

  T res;
  if (is_double<T>::value) {
    res = cblas_dasum(n, (const double *)x.data(), incx);
  } else if (is_float<T>::value) {
    res = cblas_sasum(n, (const float *)x.data(), incx);
  } else if (is_complex_double<T>::value) {
    res = cblas_dzasum(n, reinterpret_cast<const double *>(x.data()), incx);
  } else if (is_complex_float<T>::value) {
    res = cblas_scasum(n, reinterpret_cast<const float *>(x.data()), incx);
  } else {
    err_quit("blas_asum(): unsupported element type.");
  }

  return res;
}

/// @}
/// @} BLAS INTERFACE

}  // namespace slab

#endif
