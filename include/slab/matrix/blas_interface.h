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

#include "slab/matrix/config.h"
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
#include "slab/matrix/blas/dotu.h"
#include "slab/matrix/blas/iamax.h"
#include "slab/matrix/blas/nrm2.h"
#include "slab/matrix/blas/scal.h"
#include "slab/matrix/blas/sdot.h"
#include "slab/matrix/blas/swap.h"

// BLAS Level 2 Routines and Functions
#include "slab/matrix/blas/gemv.h"
#include "slab/matrix/blas/spr.h"
#include "slab/matrix/blas/spr2.h"

// BLAS Level 3 Routines and Functions
#include "slab/matrix/blas/gemm.h"

namespace slab {

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

}  // namespace slab

#endif  // SLAB_MATRIX_BLAS_INTERFACE_H_
