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

/// @file nrm2.h
/// @brief Computes the Euclidean norm of a vector.

#ifndef SLAB_MATRIX_BLAS_NRM2_H_
#define SLAB_MATRIX_BLAS_NRM2_H_

namespace slab {

/// @brief Computes the Euclidean norm of a vector
template <typename T>
inline double blas_nrm2(const Matrix<T, 1> &x) {
  double res = 0.0;

  const int n = x.size();
  const int incx = x.descriptor().strides[0];

  if (is_double<T>::value) {
    res =
        cblas_dnrm2(n, (const double *)(x.data() + x.descriptor().start), incx);
  } else if (is_float<T>::value) {
    res =
        cblas_snrm2(n, (const float *)(x.data() + x.descriptor().start), incx);
  } else if (is_complex_double<T>::value) {
    res = cblas_dznrm2(
        n, reinterpret_cast<const double *>(x.data() + x.descriptor().start),
        incx);
  } else if (is_complex_float<T>::value) {
    res = cblas_scnrm2(
        n, reinterpret_cast<const float *>(x.data() + x.descriptor().start),
        incx);
  } else {
    err_quit("blas_nrm2(): unsupported element type.");
  }

  return res;
}

}  // namespace slab

#endif
