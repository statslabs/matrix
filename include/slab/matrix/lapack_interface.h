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

/// @file lapack_interface.h
/// @brief LAPACK interface

#ifndef _SLAB_MATRIX_LAPACK_INTERFACE_H
#define _SLAB_MATRIX_LAPACK_INTERFACE_H

#include <cassert>
#include <cstddef>

#include <complex>
#include <memory>  // std::addressof

#include "slab/__config"
#include "slab/__error"

#if defined(_SLAB_USE_NO_LAPACK)
#elif defined(_SLAB_USE_MKL)
#include "mkl.h"
#else
extern "C" {
#include "lapacke.h"
}
#endif

#include "slab/matrix/matrix.h"
#include "slab/matrix/matrix_base.h"
#include "slab/matrix/traits.h"

// LAPACK Linear Equation Computational Routines
#include "slab/matrix/lapack/getrf.h"
#include "slab/matrix/lapack/potrf.h"

// LAPACK Linear Equation Driver Routines
#include "slab/matrix/lapack/gesv.h"

// Linear Least Squares (LLS) Problems: LAPACK Driver Routines
#include "slab/matrix/lapack/gels.h"

// Singular Value Decomposition: LAPACK Driver Routines
#include "slab/matrix/lapack/gesvd.h"

#endif  // _SLAB_MATRIX_LAPACK_INTERFACE_H
