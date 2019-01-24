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

/// @file copy.h
/// @brief Copies vector to another vector.

#ifndef SLAB_MATRIX_BLAS_COPY_H_
#define SLAB_MATRIX_BLAS_COPY_H_

namespace slab {

/// @brief Copies vector to another vector
template <typename T>
inline void blas_copy(const Matrix<T, 1> &x, Matrix<T, 1> &y) {
  y.clear();
  y = Matrix<T, 1>(x.size());

  const int incx = x.descriptor().strides[0];
  const int incy = y.descriptor().strides[0];

  if (is_double<T>::value)
    cblas_dcopy(x.size(), (const double *)(x.data() + x.descriptor().start),
                incx, (double *)(y.data() + y.descriptor().start), incy);
  else if (is_float<T>::value) {
    cblas_scopy(x.size(), (const float *)(x.data() + x.descriptor().start),
                incx, (float *)(y.data() + y.descriptor().start), incy);
  } else if (is_complex_double<T>::value) {
    cblas_zcopy(x.size(),
                reinterpret_cast<const double *>(x.data() + x.descriptor().start),
                incx, reinterpret_cast<double *>(y.data() + y.descriptor().start),
                incy);
  } else if (is_complex_float<T>::value) {
    cblas_ccopy(x.size(),
                reinterpret_cast<const float *>(x.data() + x.descriptor().start),
                incx, reinterpret_cast<float *>(y.data() + y.descriptor().start),
                incy);
  } else {
    err_quit("blas_copy(): unsupported element type.");
  }
}
 
}  // namespace slab

#endif
