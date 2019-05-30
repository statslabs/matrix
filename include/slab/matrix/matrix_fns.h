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

/// @file matrix_fns.h
/// @brief A collection of functions for Matrix

#ifndef _SLAB_MATRIX_MATRIX_FNS_H
#define _SLAB_MATRIX_MATRIX_FNS_H

#include <algorithm>
#include <complex>
#include <type_traits>

#include "slab/__config"
#include "slab/matrix/lapack_interface.h"

#include "slab/matrix/fns/eye.h"
#include "slab/matrix/fns/ones.h"
#include "slab/matrix/fns/zeros.h"

#include "slab/matrix/fns/as_scalar.h"
#include "slab/matrix/fns/diag.h"
#include "slab/matrix/fns/join.h"
#include "slab/matrix/fns/kron.h"
#include "slab/matrix/fns/prod.h"
#include "slab/matrix/fns/reshape.h"
#include "slab/matrix/fns/sum.h"
#include "slab/matrix/fns/trans.h"
#include "slab/matrix/fns/vectorise.h"

#include "slab/matrix/fns/misc.h"
#include "slab/matrix/fns/trig.h"

#include "slab/matrix/fns/chol.h"
#include "slab/matrix/fns/inv.h"
#include "slab/matrix/fns/lu.h"
#include "slab/matrix/fns/pinv.h"
#include "slab/matrix/fns/solve.h"

#endif  // _SLAB_MATRIX_MATRIX_FNS_H
