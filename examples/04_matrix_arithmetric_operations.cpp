#include "slab/matrix.h"
using namespace slab;

#include <iostream>
using namespace std;

int main() {
  Matrix<int, 2> m1{{1, 2, 3}, {4, 5, 6}};
  Matrix<int, 2> m2{m1};
  m1 *= 2;
  Matrix<int, 2> m3 = m1 + m2;
  Matrix<int, 2> m4{{1, 2}, {3, 4}, {5, 6}};
  Matrix<int, 2> v = matmul(m1, m4);

  cout << "\nm1 = " << m1
       << "\nm2 = " << m2
       << "\nm3 = " << m3
       << "\nm4 = " << m4
       << "\nv  = " << v
       << endl;

  return 0;
}
