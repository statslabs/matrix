# StatsLabs.Matrix - The C++ Linear Algebra Library #

The repository contains the StatsLabs.Matrix C++ linear algebra library code. StatsLabs.Matrix is the fundamental package of StatsLabs for statistical computing in C++.
 
## About Matrix

StatsLabs.Matrix is the fundamental package of StatsLabs for statistical computing in C++. The StatsLabs.Matrix library code is based on the matrix design chapter in 'The C++ Programming Language (4th Edition)' and provides:
  + A Matrix Template: Construction and Assignment; Subscripting and Slicing
  + Matrix arithmetic operations: Scalar Operations; Additions; Multiplication
  + Matrix Implementation: slice; MatrixSlice; MatrixRef; Matrix List Initialization; Matrix Access; Zero-Dimensional Matrix
  + An interface to Intel(R) MKL BLAS operations which apply to the Matrix template
  
## Example Program
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