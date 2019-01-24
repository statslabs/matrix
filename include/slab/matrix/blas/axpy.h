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

/// @file axpy.h
/// @brief Computes a vector-scalar product and adds the result to a vector.

#ifndef SLAB_MATRIX_BLAS_AXPY_H_
#define SLAB_MATRIX_BLAS_AXPY_H_

namespace slab {

/// @brief Computes a vector-scalar product and adds the result to a vector
template <typename T>
inline void blas_axpy(const T &a, const MatrixBase<T, 1> &x,
                      MatrixBase<T, 1> &y) {
  assert(x.size() == y.size());

  const int n = x.size();
  const int incx = x.descriptor().strides[0];
  const int incy = y.descriptor().strides[0];

  if (is_double<T>::value) {
    cblas_daxpy(n, a, (const double *)(x.data() + x.descriptor().start), incx,
                (double *)(y.data() + y.descriptor().start), incy);
  } else if (is_float<T>::value) {
    cblas_saxpy(n, a, (const float *)(x.data() + x.descriptor().start), incx,
                (float *)(y.data() + y.descriptor().start), incy);
  } else if (is_complex_double<T>::value) {
    cblas_zaxpy(
	n, reinterpret_cast<const double *>(&a), reinterpret_cast<const double *>(x.data() + x.descriptor().start),
        incx, reinterpret_cast<double *>(y.data() + y.descriptor().start), incy);
  } else if (is_complex_float<T>::value) {
    cblas_caxpy(
	n, reinterpret_cast<const float *>(&a), reinterpret_cast<const float *>(x.data() + x.descriptor().start),
        incx, reinterpret_cast<float *>(y.data() + y.descriptor().start), incy);
  } else {
    err_quit("blas_axpy(): unsupported element type.");
  }
}

}  // namespace slab

#endif
