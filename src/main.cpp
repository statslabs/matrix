//
// Created by Yi Pan (Institute of Cancer and Genomic Sciences) on 06/03/2018.
//

#include <cstdio>
#include <iostream>
#include "slab/matrix.h"

using namespace std;

int main() {
  slab::vec x1 = {1, 2, 3};
  slab::mat x2 = {{1, 2, 3}};
  cout << slab::matmul(x1, x2) << endl;

  slab::cube x3 = slab::ones<slab::cube>(3,3,3);
  cout << "x3 = " << x3 << endl;

  slab::mat x4 = slab::eye<slab::mat>(3,3);
  cout << "x4 = " << x4 << endl;

  return 0;
}