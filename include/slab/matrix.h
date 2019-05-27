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

/// @file matrix.h
/// @brief main header file for Statslabs::Matrix

#ifndef _SLAB_MATRIX_H
#define _SLAB_MATRIX_H

#include "slab/__config"
#include "slab/__error"

#include "slab/matrix/slice.h"

#ifdef _SLAB_USE_RCPP_AS_WRAP
#include <RcppCommon.h>
#endif

#include "slab/matrix/matrix.h"
#include "slab/matrix/matrix_ops.h"
#include "slab/matrix/packed_matrix.h"

#include "slab/matrix/matrix_fns.h"
#include "slab/matrix/type_alias.h"

#ifndef _SLAB_USE_NO_BLAS
#include "slab/matrix/blas_interface.h"
#endif

#ifndef _SLAB_USE_NO_LAPACK
#include "slab/matrix/lapack_interface.h"  // TODO: remove this header, provide BLAS interface only
#endif

#ifdef _SLAB_USE_RCPP_AS_WRAP
#include <Rcpp.h>
#endif

#endif  // _SLAB_MATRIX_H
