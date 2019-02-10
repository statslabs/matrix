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

/// @file dotc.h
/// @brief C++ template wrapper for C functions cblas_?dotc

#ifndef SLAB_MATRIX_BLAS_DOTC_H_
#define SLAB_MATRIX_BLAS_DOTC_H_

namespace slab {

/// @addtogroup blas_interface BLAS Interface
/// @{

/// @addtogroup blas_level1 BLAS Level 1
/// @{
  
/// @brief Computes a dot product of a conjugated vector with another vector.
template <typename T>
inline void blas_dotc_sub(const Matrix<T, 1> &x, const Matrix<T, 1> &y, Matrix<T, 1> &dotc) {
  assert(x.size() == y.size());

  const int n = x.size();
  const int incx = x.descriptor().strides[0];
  const int incy = y.descriptor().strides[0];
  
  dotc = Matrix<T, 1>(n);
  if (is_complex_double<T>::value) {
    cblas_zdotc_sub(n, reinterpret_cast<const double *>(x.data() + x.descriptor().start), incx,
		    reinterpret_cast<const double *>(y.data() + y.descriptor().start), incy,
		    reinterpret_cast<double *>(dotc.data()));
  } else if (is_complex_float<T>::value) {
        cblas_zdotc_sub(n, reinterpret_cast<const float *>(x.data() + x.descriptor().start), incx,
		    reinterpret_cast<const float *>(y.data() + y.descriptor().start), incy,
		    reinterpret_cast<float *>(dotc.data()));
  }
}

/// @} BLAS Level 1
/// @} BLAS Interface
 
}  // namespace slab

#endif
