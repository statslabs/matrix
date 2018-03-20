//
// Created by Yi Pan (Institute of Cancer and Genomic Sciences) on 06/03/2018.
//

#include <cstdio>
#include <iostream>
#include "slab/matrix.h"

using namespace std;

int main() {
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

  slab::vec c = {1.0, 2.0, 3.0};
  slab::vec cc = {1.0, 2.0, 3.0};
  cout << "dot(c, cc) = " << slab::blas_dot(c, cc) << endl;

  slab::fvec c2 = {4.0, 5.0, 6.0};

  slab::blas_scal(0.1, c);
  slab::blas_scal(0.1f, c2);

  cout << "c = " << c << endl;
  cout << "c2 = " << c2 << endl;

  cout << "nrm2(c) = " << slab::blas_nrm2(c) << endl;
  cout << "nrm2(c2) = " << slab::blas_nrm2(c2) << endl;

  cout << "asum(c) = " << slab::blas_asum(c) << endl;
  cout << "asum(c2) = " << slab::blas_asum(c2)  << endl;
  slab::blas_axpy(2.0, c, c);
  cout << "after axpy(2.0, c, c), c = " <<  c << endl;
  slab::blas_axpy(2.0f, c2, c2);
  cout << "after axpy(2.0, c2, c2), c2 = " <<  c2 << endl;

  slab::mat tmp = {
      {-0.63, 1.60, 0.49},
      {0.18, 0.33, 0.74},
      {-0.84, -0.82, 0.58}
  };

  //slab::ivec tmp2(3);
  slab::ivec tmp2;

  slab::lapack_getrf(tmp, tmp2);
  cout << "tmp = " << tmp << endl;
  cout << "tmp2 = " << tmp2 << endl;

  cout << "min(c) = " << slab::blas_iamax(c) << endl;
  cout << "min(c2) = " << slab::blas_iamax(c2) << endl;


  return 0;
}