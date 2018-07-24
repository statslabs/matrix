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

  return 0;
}