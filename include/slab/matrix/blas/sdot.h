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

/// @file sdot.h
/// @brief Computes a vector-vector dot product with double precision.

#ifndef SLAB_MATRIX_BLAS_SDOT_H_
#define SLAB_MATRIX_BLAS_SDOT_H_

namespace slab {

/// @brief Computes a vector-vector dot product with double precision
inline float blas_sdsdot(const float sb, const Matrix<float, 1> &sx,
                         const Matrix<float, 1> &sy) {
  assert(sx.size() == sy.size());

  const int n = sx.size();
  const int incx = sx.descriptor().strides[0];
  const int incy = sy.descriptor().strides[0];

  return cblas_sdsdot(n, sb, sx.data() + sx.descriptor().start, incx,
                      sy.data() + sy.descriptor().start, incy);
}

/// @brief Computes a vector-vector dot product with double precision
inline double blas_dsdot(const Matrix<float, 1> &sx,
                         const Matrix<float, 1> &sy) {
  assert(sx.size() == sy.size());

  const int n = sx.size();
  const int incx = sx.descriptor().strides[0];
  const int incy = sy.descriptor().strides[0];

  return cblas_dsdot(n, sx.data() + sx.descriptor().start, incx,
                     sy.data() + sy.descriptor().start, incy);
}

}  // namespace slab

#endif
