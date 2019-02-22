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

/// @file reshape.h
/// @brief change size while keeping elements.

#ifndef SLAB_MATRIX_FNS_RESHAPE_H_
#define SLAB_MATRIX_FNS_RESHAPE_H_

namespace slab {

template <typename T, std::size_t N, typename... Args>
inline auto reshape(const Matrix<T, N> &x, Args... args)
    -> decltype(Matrix<T, sizeof...(args)>()) {
  Matrix<T, sizeof...(args)> res(args...);
  std::copy(x.begin(), x.end(), res.begin());

  return res;
}

}  // namespace slab

#endif  // SLAB_MATRIX_FNS_RESHAPE_H_
