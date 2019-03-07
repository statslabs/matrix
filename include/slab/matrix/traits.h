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
// traits.h
// -----------------------------------------------------------------------------
//
#ifndef SLAB_MATRIX_TRAITS_H_
#define SLAB_MATRIX_TRAITS_H_

#include <complex>
#include <type_traits>  // std::enable_if/is_convertible

namespace slab {

template <bool B, typename T = void>
using Enable_if = typename std::enable_if<B, T>::type;

template <typename X, typename Y>
constexpr bool Same() {
  return std::is_same<X, Y>::value;
}

template <typename X, typename Y>
constexpr bool Convertible() {
  return std::is_convertible<X, Y>::value;
}

constexpr bool All() { return true; }

template <typename... Args>
constexpr bool All(bool b, Args... args) {
  return b && All(args...);
}

constexpr bool Some() { return false; }

template <typename... Args>
constexpr bool Some(bool b, Args... args) {
  return b || Some(args...);
}

struct substitution_failure {};

template <typename T>
struct substitution_succeeded : std::true_type {};

template <>
struct substitution_succeeded<substitution_failure> : std::false_type {};

template <typename M>
struct get_matrix_type_result {
  template <typename T, size_t N, typename = Enable_if<(N >= 1)>>
  static bool check(const Matrix<T, N> &m);

  template <typename T, size_t N, typename = Enable_if<(N >= 1)>>
  static bool check(const MatrixRef<T, N> &m);

  static substitution_failure check(...);

  using type = decltype(check(std::declval<M>()));
};

template <typename T>
struct has_matrix_type
    : substitution_succeeded<typename get_matrix_type_result<T>::type> {};

template <typename M>
constexpr bool Has_matrix_type() {
  return has_matrix_type<M>::value;
}

template <typename M>
using Matrix_type_result = typename get_matrix_type_result<M>::type;

template <typename M>
constexpr bool Matrix_type() {
  return Has_matrix_type<M>();
}

template <typename C>
using Value_type = typename C::value_type;

template <typename T>
struct is_double : public std::false_type {};

template <>
struct is_double<double> : public std::true_type {};

template <typename T>
struct is_float : public std::false_type {};

template <>
struct is_float<float> : public std::true_type {};

template <typename T>
struct is_complex_double : public std::false_type {};

template <>
struct is_complex_double<std::complex<double>> : public std::true_type {};

template <typename T>
struct is_complex_float : public std::false_type {};

template <>
struct is_complex_float<std::complex<float>> : public std::true_type {};

template <typename T>
void ignore(T &&) {}

}  // namespace slab

#endif  // SLAB_MATRIX_TRAITS_H_