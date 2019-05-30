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

/// @file blas_interface.h
/// @brief BLAS interface

#ifndef _SLAB_MATRIX_BLAS_INTERFACE_H
#define _SLAB_MATRIX_BLAS_INTERFACE_H

#include <cassert>
#include <cstddef>

#include <complex>
#include <memory>  // std::addressof

#include "slab/__config"
#include "slab/__error"

#if defined(_SLAB_USE_NO_BLAS)
#elif defined(_SLAB_USE_MKL)
#include "mkl.h"
#elif defined(_SLAB_USE_SUNPERF)
#include "sunperf.h"
#else
#include "cblas.h"
#endif

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
#include "slab/matrix/blas/rot.h"
#include "slab/matrix/blas/rotg.h"
#include "slab/matrix/blas/rotm.h"
#include "slab/matrix/blas/rotmg.h"
#include "slab/matrix/blas/scal.h"
#include "slab/matrix/blas/sdot.h"
#include "slab/matrix/blas/swap.h"

// BLAS Level 2 Routines and Functions
#include "slab/matrix/blas/gemv.h"
#include "slab/matrix/blas/spr.h"
#include "slab/matrix/blas/spr2.h"

// BLAS Level 3 Routines and Functions
#include "slab/matrix/blas/gemm.h"

#endif  // _SLAB_MATRIX_BLAS_INTERFACE_H
