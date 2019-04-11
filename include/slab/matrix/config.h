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

/// @file config.h
/// @brief Configuration Header

#ifndef STATSLABS_MATRIX_CONFIG_H_
#define STATSLABS_MATRIX_CONFIG_H_

#if !defined(USE_MKL)
// #define USE_MKL
//// Uncomment the above line if using Intel Math Kernel Library (Intel MKL)
#endif

#if !defined(USE_SUNPERF)
// #define USE_SUNPERF
//// Uncomment the above line if using Sun Performance Library
#endif

#if !defined(USE_R_BLAS)
// #define USE_R_BLAS
//// Uncomment the above line if using R's internal BLAS
#endif

#if !defined(USE_R_LAPACK)
// #define USE_R_LAPACK
//// Uncomment the above line if using R's internal LAPACK
#endif

#if !defined(USE_R_ERROR)
// #define USE_R_ERROR
//// Uncomment the above line if using Rf_error/Rf_warning in file R_Ext/Error.h
#endif

#if !defined(USE_RCPP_AS_WRAP)
// #define USE_RCPP_AS_WRAP
//// Uncomment the above line if using Rcpp::as (for conversion of objects from
//// R to C++) and Rcpp::wrap (for conversion from C++ to R)
#endif

#if defined(USE_MKL)
#include "mkl.h"
#elif defined(USE_SUNPERF)
#include "sunperf.h"
#else
extern "C" {
#include "cblas.h"
}
#endif

#if defined (USE_R_ERROR)
#  include <R_ext/Error.h>
#  define SLAB_ERROR(...) Rf_error(__VA_ARGS__)
#  define SLAB_WARNING(...) Rf_warning(__VA_ARGS__)
#else
#  include "slab/matrix/error.h"
#  define SLAB_ERROR(...) err_quit(__VA_ARGS__)
#  define SLAB_WARNING(...) err_msg(__VA_ARGS__)
#endif

#endif  // STATSLABS_MATRIX_CONFIG_H_
