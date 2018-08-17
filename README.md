# Statslabs.Matrix - The C++ Linear Algebra Library #

The repository contains the Statslabs.Matrix C++ linear algebra library code. Statslabs.Matrix is the fundamental package of Statslabs for statistical computing in C++.

## About Matrix

Statslabs.Matrix is the fundamental package of Statslabs for statistical computing in C++. The Statslabs.Matrix library code is based on the matrix design chapter in 'The C++ Programming Language (4th Edition)' and provides:
  + A Matrix Template: Construction and Assignment; Subscripting and Slicing
  + Matrix arithmetic operations: Scalar Operations; Additions; Multiplication
  + Matrix Implementation: slice; MatrixSlice; MatrixRef; Matrix List Initialization; Matrix Access; Zero-Dimensional Matrix
  + An interface to Intel(R) MKL BLAS operations which apply to the Matrix template

## Prerequisites

    CMake >= 3.0
    Intel Math Kernel Library (Intel MKL)
   
## Installation on Ubuntu / macOS
1. Clone the repository.
   ```sh
   git clone git@github.com:statslabs/matrix.git
   ```
2. Configure the project.
   ```sh
   cd matrix
   mkdir build && cd build
   cmake ..
   ```
3. Compile and install the library.
   ```sh
   make
   make install
   ```

## Example program
`examples/demo.cc`:
```c
#include <iostream>
#include "slab/matrix.h"

using namespace std;
using namespace slab;

int main() {
  mat A = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9}
  };

  mat B = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9}
  };

  // Element-wise addition
  cout << "A + B = " << A + B << endl;

  // Element-wise subtraction
  cout << "A - B = " << A - B << endl;

  // Element-wise multiplication
  cout << "A * B = " << A * B << endl;

  // Element-wise division
  cout << "A / B = " << A / B << endl;

  // Matrix multiplication
  cout << "matmul(A, B) = " << matmul(A, B) << endl;

  return 0;
}
```

## Integration of Matrix in your own project
To make the project simple enough, we will create a CMake project for `demo.cc`.

1. Make a project folder.
   ```sh
   mkdir example && cd example
   ```

2. Create `demo.cc` and `CMakeLists.txt` in the project folder where file `CMakeLists.txt` should look like:
   ```cmake
   cmake_minimum_required(VERSION 3.0)
   project(example)
   set(CMAKE_CXX_STANDARD 11)
   add_executable(example demo.cc)

   find_package(Matrix REQUIRED)
   target_link_libraries(example Statslabs::Matrix)
   ```

3. Perform a out-of-source build.
   ```sh
   mkdir build && cd build
   cmake ..
   make
   ```

4. Run the program.
   ```sh
   ./example
   ```