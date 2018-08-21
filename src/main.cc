//
// Created by Yi Pan (Institute of Cancer and Genomic Sciences) on 06/03/2018.
//

#include <complex>
#include <iostream>
#include "slab/matrix.h"

using namespace std;

int main() {
  slab::mat a = {{6.80, -6.05, -0.45, 8.32, -9.67},
                 {-2.11, -3.30, 2.58, 2.71, -5.14},
                 {5.66, 5.36, -2.70, 4.35, -7.26},
                 {5.97, -4.44, 0.27, -7.17, 6.08},
                 {8.23, 1.08, 9.04, 2.14, -6.87}};

  slab::mat b = {{4.02, -1.56, 9.81},
                 {6.19, 4.00, -4.09},
                 {-8.22, -8.67, -4.57},
                 {-7.57, 1.75, -8.61},
                 {-3.03, 2.86, 8.99}};

  slab::ivec ipiv(5);

  int info = slab::lapack_gesv(a, ipiv, b);

  cout << "a = " << a << endl;
  cout << "b = " << b << endl;

  slab::mat A = {{1, 2, 3}, {0, 1, 4}, {5, 6, 0}};

  cout << "inverse of A is " << solve(A) << endl;

  slab::mat B;
  slab::inv(B, A);

  cout << "inverse of A is " << B << endl;

  slab::mat x = {{1, 2, 3, 4}, {5, 6, 7, 8}, {3, 4, 5, 6}};

  slab::mat x_inv;

  slab::pinv(x_inv, x);
  cout << "x_inv = " << x_inv << endl;

  return 0;
}
