#include <iostream>
using namespace std;

#include "slab/matrix.h"
using namespace slab;

int main() {

  Matrix<int, 2> m{
      {01, 02, 03},
      {11, 12, 13}
  };

  cout << "m = " << m << endl;

  m(1, 2) = 99;         // overwrite the element in row 1 column 2; that is 13
//  auto d1 = m(1);      // error: too few subscripts
//  auto d2 = m(1,2,3);  // error: too many subscripts

  cout << "m = " << m << endl;

  Matrix<int, 2> m2{
      {01, 02, 03},
      {11, 12, 13},
      {21, 22, 23}
  };

  auto m22 = m2(slice{1, 2}, slice{0, 3});

  cout << "m2 = " << m2 << endl;
  cout << "m22 = " << m22 << endl;

  m2(slice{1, 2}, slice{0, 3}) = {
      {111, 112, 113},
      {121, 122, 123}
  };

  cout << "m2 = " << m2 << endl;

  Matrix<int, 2> m3 = {
      {01, 02, 03},
      {11, 12, 13},
      {21, 22, 23}
  };

  auto m31 = m3(slice{1, 2}, 1);  // m31 becomes {{12},{22}}
  auto m32 = m3(slice{1, 2}, 0);  // m32 becomes {{11},{21}}
  auto x = m3(1, 2);              // x == 13

  cout << "m31 = " << m31 << endl;
  cout << "m32 = " << m32 << endl;
  cout << "x = " << x << endl;

  return 0;
}
