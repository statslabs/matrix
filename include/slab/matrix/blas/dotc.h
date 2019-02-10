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

/// @file dot.h
/// @brief Computes a vector-vector dot product.

#ifndef SLAB_MATRIX_BLAS_DOT_H_
#define SLAB_MATRIX_BLAS_DOT_H_

namespace slab {

/// @brief Computes a vector-vector dot product
template <typename T>
inline T blas_dot(const Matrix<T, 1> &x, const Matrix<T, 1> &y) {
  assert(x.size() == y.size());

  const int n = x.size();
  const int incx = x.descriptor().strides[0];
  const int incy = y.descriptor().strides[0];

  T res = 0.0;
  if (is_double<T>::value) {
    res = cblas_ddot(n, (const double *)(x.data() + x.descriptor().start), incx,
                     (const double *)(y.data() + y.descriptor().start), incy);
  } else if (is_float<T>::value) {
    res = cblas_sdot(n, (const float *)(x.data() + x.descriptor().start), incx,
                     (const float *)(y.data() + y.descriptor().start), incy);
  } else {
    err_quit("blas_dot(): unsupported element type.");
  }

  return res;
}

}  // namespace slab

#endif
