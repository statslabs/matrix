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

/// @file ones.h
/// @brief generate object filled with ones.

#ifndef STATSLABS_MATRIX_ONES_H_
#define STATSLABS_MATRIX_ONES_H_

namespace slab {

template <typename M, typename... Args>
Enable_if<Matrix_type<M>(), M> ones(Args... args) {
  assert(M::order() == sizeof...(args));
  M res(args...);
  res = 1;

  return res;
}

}  // namespace slab

#endif  // STATSLABS_MATRIX_ONES_H_
