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
/// @brief C++ wrapper for C functions cblas_?sdot

#ifndef SLAB_MATRIX_BLAS_SDOT_H_
#define SLAB_MATRIX_BLAS_SDOT_H_

namespace slab {

/// @addtogroup blas_interface BLAS Interface
/// @{

/// @addtogroup blas_level1 BLAS Level 1
/// @{

/// @brief Computes a vector-vector dot product with double precision.
///
/// The sdsdot routine computes the inner product of two vectors with
/// double precision, adds scalar value \f$sb\f$ to the inner product,
/// and outputs the final results in single precision.
///
/// @param sb Scalar with type double.
/// @param sx Vector with type fvec.
/// @param sy Vector with type fvec.
/// @return The result of the dot product of sx and sy (with sb added).
///
inline float blas_sdsdot(const float sb, const Matrix<float, 1> &sx,
                         const Matrix<float, 1> &sy) {
  assert(sx.size() == sy.size());

  const std::size_t n = sx.size();
  const std::size_t incx = sx.descriptor().strides[0];
  const std::size_t incy = sy.descriptor().strides[0];

  float res = cblas_sdsdot(
      (const int)n, (const float)sb,
      (const float *)(sx.data() + sx.descriptor().start), (const int)incx,
      (const float *)(sy.data() + sy.descriptor().start), (const int)incy);

  return res;
}

/// @brief Computes a vector-vector dot product with double precision.
///
/// The dsdot routine computes the inner product of two vectors with
/// double precision, and outputs the final results in double
/// precision.
///
/// @param sx Vector with type fvec
/// @param sy Vector with type fvec
/// @return The result of the dot product of sx and sy
///
inline double blas_dsdot(const Matrix<float, 1> &sx,
                         const Matrix<float, 1> &sy) {
  assert(sx.size() == sy.size());

  const std::size_t n = sx.size();
  const std::size_t incx = sx.descriptor().strides[0];
  const std::size_t incy = sy.descriptor().strides[0];

  double res = cblas_dsdot(
      (const int)n, (const float *)(sx.data() + sx.descriptor().start),
      (const int)incx, (const float *)(sy.data() + sy.descriptor().start),
      (const int)incy);

  return res;
}

/// @} BLAS Level 1
/// @} BLAS Interface

}  // namespace slab

#endif
