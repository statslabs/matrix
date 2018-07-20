/*******************************************************************************
 * EXPECTED OUTPUT:
 *******************************************************************************

d1 = 1
d2 = 2
e1 = 100
e2 = 50
e2a = 6000
s1 = 100
s2 = 300000

m = {{0,1,2,3},{10,11,12,13},{20,21,22,23}}

v = {10,11,12,13}

dval1 = 12
dval2 = 12
dval3 = 12

 ******************************************************************************/

#include <iostream>
using namespace std;

#include "slab/matrix.h"
using namespace slab;

int main() {
  Matrix<double, 0> m0{1};           // zero dimensions: a scalar
  Matrix<double, 1> m1{1, 2, 3, 4};  // one dimension: a vector (4 elements)
  Matrix<double, 2> m2{              // two dimensions (4*3 elements)
      {00, 01, 02, 03},  // row 0
      {10, 11, 12, 13},  // row 1
      {20, 21, 22, 23}   // row 2
  };

  // three dimensions (4*7*9 elements), all 0-initialized
  Matrix<double, 3> m3(4, 7, 9);
  // 17 dimensions (no elements so far)
  Matrix<complex<double>, 17> m17;

  Matrix<double, 2> md;  // OK
  Matrix<string, 2> ms;  // OK: just don't try arithmetic operations

  // 3-by-2 matrix of 2-by-2 matrices
  // a matrix is a plausible "number"
  Matrix<Matrix<int, 2>, 2> mm{
      {  // row 0
          {{1, 2}, {3, 4}},  // col0
          {{4, 5}, {6, 7}},  // col1
      },
      {  // row 1
          {{8, 9}, {0, 1}},  // col0
          {{2, 3}, {4, 5}},  // col1
      },
      {  // row 2
          {{1, 2}, {3, 4}},  // col0
          {{4, 5}, {6, 7}},  // col1
      }
  };

//  Matrix<char, 2> mc1(2, 3, 4);  // error: two many dimension sizes

//  Matrix<char, 2> mc2{
//      {'1', '2', '3'}
//  };

//  Matrix<char, 2> mc3{
//      {'1', '2', '3'},
//      {'4', '5'}  // error: element missing for third column
//  };

  Matrix<double, 1> m1_new(100);       // one dimension: a vector (100 elements)
  Matrix<double, 2> m2_new(50, 6000);  // two dimensions: 50*6000

  auto d1 = m1_new.order();            // 1
  auto d2 = m2_new.order();            // 2

  auto e1 = m1_new.extent(0);          // 100
//  auto e1a = m1_new.extent(1);         // error: m1 is one-dimensional

  auto e2 = m2_new.extent(0);          // 50
  auto e2a = m2_new.extent(1);         // 6000

  auto s1 = m1_new.size();             // 100
  auto s2 = m2_new.size();             // 50*6000

  cout << "\nd1 = " << d1 << "\nd2 = " << d2
       << "\ne1 = " << e1
       << "\ne2 = " << e2 << "\ne2a = " << e2a
       << "\ns1 = " << s1 << "\ns2 = " << s2
       << endl;

  Matrix<double, 2> m{                 // two dimensions (4*3 elements)
      {00, 01, 02, 03},  // row 0
      {10, 11, 12, 13},  // row 1
      {20, 21, 22, 23}   // row 2
  };

  cout << "\nm = " << m << endl;

  Matrix<double, 1> v = m[1];
  cout << "\nv = " << v << endl;

  double dval1 = m(1, 2);              // 12
  double dval2 = m[1][2];              // 12
  double dval3 = v[2];                 // 12

  cout << "\ndval1 = " << dval1
       << "\ndval2 = " << dval2
       << "\ndval3 = " << dval3
       << endl;

  return 0;
}
