//
// Copyright 2018 The Statslabs Authors.
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

/// @file blas_interface.h
/// @brief BLAS interface

#ifndef SLAB_MATRIX_BLAS_INTERFACE_H_
#define SLAB_MATRIX_BLAS_INTERFACE_H_

#include <cassert>
#include <cstddef>

#include <complex>

#ifdef USE_MKL
#include "mkl.h"
#else
extern "C" {
#include "cblas.h"
}
#endif

#include "slab/matrix/error.h"
#include "slab/matrix/matrix.h"
#include "slab/matrix/matrix_base.h"
#include "slab/matrix/traits.h"

// BLAS Level 1 Routines and Functions
#include "slab/matrix/blas/asum.h"
#include "slab/matrix/blas/axpy.h"
#include "slab/matrix/blas/copy.h"
#include "slab/matrix/blas/dot.h"
#include "slab/matrix/blas/dotc.h"
#include "slab/matrix/blas/nrm2.h"
#include "slab/matrix/blas/sdot.h"

// BLAS Level 2 Routines and Functions
#include "slab/matrix/blas/gemv.h"

// BLAS Level 3 Routines and Functions
#include "slab/matrix/blas/gemm.h"

namespace slab {

/// @addtogroup blas_interface BLAS INTERFACE
/// @{

/// @addtogroup blas_level1 BLAS Level 1
/// @{

// Computes the parameters for a Givens rotation
// template<typename T>
// void blas_rotg(Matrix<T, 1> &a, Matrix<T, 1> &b, Matrix<T, 1> &c, Matrix<T,
// 1> &s) {
//  if (is_double<T>::value) {
//    cblas_drotg((double *) a.data(),
//                (double *) b.data(),
//                (double *) c.data(),
//                (double *) s.data());
//  } else if (is_float<T>::value) {
//    cblas_srotg((float *) a.data(), (float *) b.data(), (float *) c.data(),
//    (float *) s.data());
//  } else if (is_complex_double<T>::value) {
//    cblas_zrotg((std::complex<double> *) a.data(), (const std::complex<double>
//    *) b.data(),
//                (double *) c.data(), (std::complex<double> *) s.data());
//  } else if (is_complex_float<T>::value) {
//    cblas_crotg((std::complex<float> *) a.data(), (const std::complex<float>
//    *) b.data(),
//                (float *) c.data(), (std::complex<float> *) s.data());
//  }
//}

/// @brief Computes the product of a vector by a scalar
template <typename T>
inline void blas_scal(const T a, Matrix<T, 1> &x) {
  const int n = x.size();
  const int incx = x.descriptor().strides[0];

  if (is_double<T>::value) {
    cblas_dscal(n, (const double)a, (double *)(x.data() + x.descriptor().start),
                incx);
  } else if (is_float<T>::value) {
    cblas_sscal(n, (const float)a, (float *)(x.data() + x.descriptor().start),
                incx);
  } else if (is_complex_double<T>::value) {
    cblas_zdscal(n, (const double)a,
                 reinterpret_cast<double *>(x.data() + x.descriptor().start),
                 incx);
  } else if (is_complex_float<T>::value) {
    cblas_csscal(n, (const float)a,
                 reinterpret_cast<float *>(x.data() + x.descriptor().start),
                 incx);
  } else {
    err_quit("blas_scal(): unsupported element type.");
  }
}

template <typename T>
inline void blas_scal(const std::complex<T> &a, Matrix<std::complex<T>, 1> &x) {
  const int n = x.size();
  const int incx = x.descriptor().strides[0];

  if (is_double<T>::value) {
    cblas_zscal(n, reinterpret_cast<const double *>(&a),
                reinterpret_cast<double *>(x.data() + x.descriptor().start),
                incx);
  } else if (is_float<T>::value) {
    cblas_cscal(n, reinterpret_cast<const float *>(&a),
                reinterpret_cast<float *>(x.data() + x.descriptor().start),
                incx);
  } else {
    err_quit("blas_scal(): unsupported element type.");
  }
}

/// @brief Swaps a vector with another vector.
///
/// @param x a vector.
/// @param y another vector.
template <typename T>
inline void blas_swap(Matrix<T, 1> &x, Matrix<T, 1> &y) {
  assert(x.size() == y.size());

  const int n = x.size();
  const int incx = x.descriptor().strides[0];
  const int incy = y.descriptor().strides[0];

  if (is_double<T>::value) {
    cblas_dswap(n, (double *)(x.data() + x.descriptor().start), incx,
                (double *)(y.data() + y.descriptor().start), incy);
  } else if (is_float<T>::value) {
    cblas_sswap(n, (float *)(x.data() + x.descriptor().start), incx,
                (float *)(y.data() + y.descriptor().start), incy);
  } else if (is_complex_double<T>::value) {
    cblas_zswap(
        n, reinterpret_cast<double *>(x.data() + x.descriptor().start), incx,
        reinterpret_cast<double *>(y.data() + y.descriptor().start), incy);
  } else if (is_complex_float<T>::value) {
    cblas_cswap(
        n, reinterpret_cast<float *>(x.data() + x.descriptor().start), incx,
        reinterpret_cast<float *>(y.data() + y.descriptor().start), incy);
  } else {
    err_quit("blas_swap(): unsupported element type.");
  }
}

/// @brief Finds the index of the element with maximum absolute value
template <typename T>
inline std::size_t blas_iamax(const Matrix<T, 1> &x) {
  std::size_t res = 0;
  std::size_t incx = x.descriptor().strides[0];

  if (is_double<T>::value) {
    res = cblas_idamax(x.size(),
                       (const double *)(x.data() + x.descriptor().start), incx);
  } else if (is_float<T>::value) {
    res = cblas_isamax(x.size(),
                       (const float *)(x.data() + x.descriptor().start), incx);
  } else if (is_complex_double<T>::value) {
    res = cblas_izamax(
        x.size(),
        reinterpret_cast<const double *>(x.data() + x.descriptor().start),
        incx);
  } else if (is_complex_float<T>::value) {
    res = cblas_icamax(
        x.size(),
        reinterpret_cast<const float *>(x.data() + x.descriptor().start), incx);
  } else {
    err_quit("blas_iamax(): unsupported element type.");
  }

  return res;
}
/// @}

/// @addtogroup blas_level2 BLAS Level 2
/// @{

template <typename T, typename TRI>
inline void blas_spr(const T &alpha, const MatrixBase<T, 1> &x,
                     SymmetricMatrix<T, TRI> &ap) {
  assert(x.size() == ap.n_rows());

  CBLAS_UPLO uplo;
  if (is_upper<TRI>::value)
    uplo = CblasUpper;
  else if (is_lower<TRI>::value)
    uplo = CblasLower;

  const int incx = x.descriptor().strides[0];

  if (is_double<T>::value) {
    cblas_dspr(CblasRowMajor, uplo, x.size(), (const double)alpha,
               (const double *)x.data(), incx, (double *)ap.data());
  } else if (is_float<T>::value) {
    cblas_sspr(CblasRowMajor, uplo, x.size(), (const float)alpha,
               (const float *)x.data(), incx, (float *)ap.data());
  } else {
    err_quit("blas_spr(): unsupported element type.");
  }
}

template <typename T, typename TRI>
inline void blas_spr2(const T &alpha, const MatrixBase<T, 1> &x,
                      const MatrixBase<T, 1> &y, SymmetricMatrix<T, TRI> &ap) {
  assert(x.size() == y.size());
  assert(x.size() == ap.n_rows());

  CBLAS_UPLO uplo;
  if (is_upper<TRI>::value)
    uplo = CblasUpper;
  else if (is_lower<TRI>::value)
    uplo = CblasLower;

  const int incx = x.descriptor().strides[0];
  const int incy = y.descriptor().strides[0];

  if (is_double<T>::value) {
    cblas_dspr2(CblasRowMajor, uplo, x.size(), (const double)alpha,
                (const double *)x.data(), incx, (const double *)y.data(), incy,
                (double *)ap.data());
  } else if (is_float<T>::value) {
    cblas_sspr2(CblasRowMajor, uplo, x.size(), (const float)alpha,
                (const float *)x.data(), incx, (const float *)y.data(), incy,
                (float *)ap.data());
  } else {
    err_quit("blas_spr2(): unsupported element type.");
  }
}

/// @}

/// @addtogroup blas_level3 BLAS Level 3
/// @{
/// @}
/// @} BLAS INTERFACE

}  // namespace slab

#endif  // SLAB_MATRIX_BLAS_INTERFACE_H_
