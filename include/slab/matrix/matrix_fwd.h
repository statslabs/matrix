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

///@file matrix_fwd.h
///@brief This file is essentially used to forward declare the main types

#ifndef SLAB_MATRIX_MATRIX_FWD_H_
#define SLAB_MATRIX_MATRIX_FWD_H_

#include <cstddef>

namespace slab {

// Declarations
struct slice;

template <std::size_t N>
struct MatrixSlice;

template <typename T, std::size_t N>
class Matrix;

template <typename T, std::size_t N>
class MatrixRef;

template <typename T>
Matrix<T, 2> transpose(const Matrix<T, 1> &a);

template <typename T>
Matrix<T, 2> transpose(const MatrixRef<T, 1> &a);

template <typename T>
Matrix<T, 2> transpose(const Matrix<T, 2> &a);

template <typename T>
Matrix<T, 2> transpose(const MatrixRef<T, 2> &a);

template <typename T>
Matrix<T, 2> inverse(const Matrix<T, 2> &a);

template <typename T>
Matrix<T, 2> inverse(const MatrixRef<T, 2> &a);

}  // namespace slab

#endif  // SLAB_MATRIX_MATRIX_FWD_H_
