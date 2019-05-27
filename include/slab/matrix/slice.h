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

/// @file slice.h
/// @brief slice struct

#ifndef _SLAB_MATRIX_SLICE_H
#define _SLAB_MATRIX_SLICE_H

#include <cstddef>
#include <iostream>

_SLAB_BEGIN_NAMESPACE

// ---------------------------------------------------------------------- //
//                                   Slice
// ---------------------------------------------------------------------- //
// A slice describes a sequence of elements in some dimension (or row) of a
// matrix. It is a triple comprised of a starting index, a number of elements,
// and the stride between subsequent elements.
//
// The special members slice::none and slice::all represent the selection of
// no or all elements in a particular dimension.
struct slice {
  slice() : start(-1), length(-1), stride(1) {}
  explicit slice(std::size_t s) : start(s), length(-1), stride(1) {}
  slice(std::size_t s, std::size_t l, std::size_t n = 1)
      : start(s), length(l), stride(n) {}

  std::size_t operator()(std::size_t i) const { return start + i * stride; }

  static slice all;

  std::size_t start;   // first index
  std::size_t length;  // number of indices included
  std::size_t stride;  // distance between elements in sequence
};

// FIXME
// slice slice::all{0, std::size_t(-1), 1};

template <typename C, typename T>
std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& os,
                                     const slice& s) {
  return os << '(' << s.start << ' ' << s.length << ' ' << s.stride << ')';
}

_SLAB_END_NAMESPACE

#endif  // _SLAB_MATRIX_SLICE_H
