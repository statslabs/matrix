//
// Copyright 2018 The StatsLabs Authors.
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
// matrix_base.h
// -----------------------------------------------------------------------------
//
#ifndef STATSLABS_MATRIX_MATRIX_BASE_H_
#define STATSLABS_MATRIX_MATRIX_BASE_H_

#include <cstddef>
#include "slab/matrix/matrix_slice.h"
#include "slab/matrix/support.h"

template<typename T, std::size_t N>
class MatrixBase {
 public:
  static constexpr std::size_t order = N;
  using value_type = T;

  MatrixBase() = default;
  MatrixBase(MatrixBase &&) = default;
  MatrixBase &operator=(MatrixBase &&) = default;
  MatrixBase(MatrixBase const &) = default;
  MatrixBase &operator=(MatrixBase const &) = default;

  template<typename... Exts>
  explicit MatrixBase(Exts... exts);

  std::size_t extent(std::size_t n) const { return desc.extents[n]; }
  std::size_t rows() const { return desc.extents[0]; }
  std::size_t cols() const { return desc.extents[1]; }
  const MatrixSlice<N> &descriptor() const { return desc; }

 private:
  MatrixSlice<N> desc;
};

template<typename T, std::size_t N>
template<typename... Exts>
MatrixBase<T, N>::MatrixBase(Exts... exts)
    : desc(exts...) {}

#endif // STATSLABS_MATRIX_MATRIX_BASE_H_
