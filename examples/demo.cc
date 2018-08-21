#include <iostream>
#include "slab/matrix.h"

using namespace std;
using namespace slab;

int main() {
  mat A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

  mat B = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

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