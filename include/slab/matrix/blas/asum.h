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

/// @file asum.h
/// @brief Computes the sum of magnitudes of the vector elements.

#ifndef SLAB_MATRIX_BLAS_ASUM_H_
#define SLAB_MATRIX_BLAS_ASUM_H_

namespace slab {

/// @brief Computes the sum of magnitudes of the vector elements
/// @param x a vector.
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

}  // namespace slab

#endif
