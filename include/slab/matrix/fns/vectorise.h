//
// Copyright 2019 The Statslabs Authors.
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

/// @file vectorise.h
/// @brief convert matrix to vector.

#ifndef SLAB_MATRIX_FNS_VECTORISE_H_
#define SLAB_MATRIX_FNS_VECTORISE_H_

namespace slab {

template <typename T, std::size_t N>
inline Matrix<T, 1> vectorise(const Matrix<T, N> &x) {
  return reshape(x, x.size());
}

}  // namespace slab

#endif  // SLAB_MATRIX_FNS_VECTORISE_H_
