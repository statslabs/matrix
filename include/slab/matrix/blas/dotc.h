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

/// @file dotc.h
/// @brief C++ template wrapper for C functions cblas_?dotc

#ifndef _SLAB_MATRIX_BLAS_DOTC_H
#define _SLAB_MATRIX_BLAS_DOTC_H

_SLAB_BEGIN_NAMESPACE

/// @addtogroup blas_interface BLAS Interface
/// @{

/// @addtogroup blas_level1 BLAS Level 1
/// @{

/// @brief Computes a dot product of a conjugated vector with another vector.
///
/// The dotc routines perform a vector-vector operation defined as:
/// \f[
/// res = \sum_{i=1}^nConjg(x_i) \times y_i
/// \f]
/// where \f$x_i\f$ and \f$y_i\f$ are elements of vectors \f$x\f$ and \f$y\f$.
///
/// @param x Vector with type cx_vec/cx_fvec.
/// @param y Vector with type cx_vec/cx_fvec.
/// @param dotc Contains the result of the dot product of the conjugated x and
/// unconjugated y.
/// @return Void.
///
template <typename T>
inline void blas_dotc_sub(const MatrixBase<T, 1> &x, const MatrixBase<T, 1> &y,
                          T &dotc) {
  _SLAB_ASSERT(x.size() == y.size(),
               "blas_dotc_sub(): incompatible vector dimensions");
  const std::size_t n = x.size();
  const std::size_t incx = x.descriptor().strides[0];
  const std::size_t incy = y.descriptor().strides[0];
  const T *x_ptr = x.data() + x.descriptor().start;
  const T *y_ptr = y.data() + y.descriptor().start;

  if (is_complex_double<T>::value) {
    cblas_zdotc_sub((const SLAB_INT)n, (const void *)x_ptr,
                    (const SLAB_INT)incx, (const void *)y_ptr,
                    (const SLAB_INT)incy, (void *)(&dotc));
  }
#ifndef _SLAB_USE_R_BLAS
  else if (is_complex_float<T>::value) {
    cblas_cdotc_sub((const SLAB_INT)n, (const void *)x_ptr,
                    (const SLAB_INT)incx, (const void *)y_ptr,
                    (const SLAB_INT)incy, (void *)(&dotc));
  }
#endif
  else {
    _SLAB_ERROR("blas_dotc_sub(): unsupported element type.");
  }
}

/// @} BLAS Level 1
/// @} BLAS Interface

_SLAB_END_NAMESPACE

#endif
