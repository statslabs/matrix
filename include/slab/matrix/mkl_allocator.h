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

#ifndef SLAB_MATRIX_MKL_ALLOCATOR_H_
#define SLAB_MATRIX_MKL_ALLOCATOR_H_

template <class T>
struct MklAllocator
{
  typedef T value_type;
  MklAllocator() noexcept {} //default ctor not required by C++ Standard Library

  // A converting copy constructor:
  template<class U> MklAllocator(const MklAllocator<U>&) noexcept {}
  template<class U> bool operator==(const MklAllocator<U>&) const noexcept
  {
    return true;
  }
  template<class U> bool operator!=(const Mallocator<U>&) const noexcept
  {
    return false;
  }
  T* allocate(const size_t n) const;
  void deallocate(T* const p, size_t) const noexcept;
};

template <class T>
T* MklAllocator<T>::allocate(const size_t n) const
{
  if (n == 0)
  {
    return nullptr;
  }
  if (n > static_cast<size_t>(-1) / sizeof(T))
  {
    throw std::bad_array_new_length();
  }
  void* const pv = malloc(n * sizeof(T));
  if (!pv) { throw std::bad_alloc(); }
  return static_cast<T*>(pv);
}

template<class T>
void MklAllocator<T>::deallocate(T * const p, size_t) const noexcept
{
  free(p);
}

#endif // SLAB_MATRIX_MKL_ALLOCATOR_H_
