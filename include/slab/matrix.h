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
// -----------------------------------------------------------------------------
// matrix.h
// -----------------------------------------------------------------------------
//
#ifndef SLAB_MATRIX_H_
#define SLAB_MATRIX_H_

#include <cassert>
#include <cerrno>   // for definition of errno
#include <cmath>
#include <cstdarg>  // ISO C variable arguments
#include <cstddef>  // std::size_t
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <array>
#include <complex>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <numeric> // std::inner_product
#include <string>
#include <type_traits> // std::enable_if/is_convertible
#include <vector>

#include "mkl.h"
// #include "slab/matrix/config.h"

namespace slab {

#include "slab/matrix/error.h"

#include "slab/matrix/matrix_fwd.h"
#include "slab/matrix/traits.h"

#include "slab/matrix/slice.h"
#include "slab/matrix/support.h"

template<typename T, std::size_t N>
using MatrixInitializer = typename matrix_impl::MatrixInit<T, N>::type;

#include "slab/matrix/matrix_slice.h"
#include "slab/matrix/matrix_ref.h"
#include "slab/matrix/matrix.h"
#include "slab/matrix/packed_matrix.h"

#include "slab/matrix/blas_interface.h"
#include "slab/matrix/lapack_interface.h"

#include "slab/matrix/matrix_ops.h"
#include "slab/matrix/type_alias.h"

} // namespace slab

#endif // SLAB_MATRIX_H_
