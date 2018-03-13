//
// Created by Yi Pan (Institute of Cancer and Genomic Sciences) on 06/03/2018.
//

#include <cstdio>
#include <iostream>
#include "slab/matrix.h"

using namespace std;

int main()
{
  slab::fmat A =
      {{1, 2, 3},
       {4, 5, 6},
       {7, 8, 9}};

  slab::fmat B =
      {{1, 2, 3},
       {4, 5, 6},
       {7, 8, 9}};

  slab::fmat C = A * B;

  cout << "C  = " << C << endl;
  cout << "AB = " << slab::matmul(A, B) << endl;
  cout << "A' = " << slab::transpose(A) << endl;

  auto D = slab::reshape(C, 1, 9);
  cout << "D = " << D << endl;

  return 0;
}