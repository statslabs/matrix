//
// Created by Yi Pan (Institute of Cancer and Genomic Sciences) on 06/03/2018.
//

#include <cstdio>
#include <complex>
#include <iostream>
#include "slab/matrix.h"

using namespace std;

int main() {
/*
  slab::vec dv = {10, 9, 8, 7, 6};

  dv.subvec(2, 4) = {1, 2, 3};

  slab::mat x1 = {
      {1, 2, 3},
      {4, 5, 6}
  };

  slab::mat x2 = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9}
  };

  slab::mat x3 = {
      {10, 10},
      {10, 10}
  };

  cout << slab::join_cols(x1, x2) << endl;
  cout << slab::join_rows(x1, x3) << endl;

  cout << x2(slab::slice{1}, slab::slice{1}) << endl;
  cout << slab::vectorise(x2(slab::slice{1}, slab::slice{1})) << endl;

  cout << dv << endl;

  cout << "x2.submat(1,1,2,2) = " << x2.submat(1,1,2,2) << endl;
  cout << "x1.submat(0,1,1,2) = " << x1.submat(0,1,1,2) << endl;

  x2.submat(1,1,2,2) = x1.submat(0,1,1,2);
  cout << "x2.submat(1,1,2,2) = " << x2.submat(1,1,2,2) << endl;

  cout << x2 << endl;
*/

//  slab::SymmetricMatrix<double, slab::upper> sm(3);
//
//  cout << "sm =\n" << sm << endl;
//
//  for (size_t i = 0; i != sm.n_rows(); ++i) {
//    for (size_t j = 0; j <= i; ++j) {
//      sm(i, j) = i * 10 + j;
//    }
//  }
//
//  cout << "sm =\n" << sm << endl;
//
//  slab::vec xvec = {1, 2, 3};
//  slab::blas_spr(2.0, xvec, sm);
//
//  cout << "sm =\n" << sm << endl;

//  slab::HermitianMatrix<std::complex<double>, slab::lower> ml(3);
//  for (unsigned i = 0; i < ml.n_rows(); ++ i) {
//    for (unsigned j = 0; j < i; ++ j)
//      ml (i, j) = std::complex<double> (3 * i + j, 3 * i + j);
//    ml (i, i) = std::complex<double> (4 * i, 0);
//  }
//
//  cout << "ml =\n" << ml << endl;


  // std::complex<double> data(2, 3);
  // cout << "da = " << data << endl;

  // std::complex<double> data2(4, 5);
  // cout << "da2 = " << data2 << endl;

  // ComplexConj<double> cc(data);
  // cout << "cc = " << std::complex<double>(cc) << endl;

  // cc = data2;

  // cout << "da = " << data << endl;


  slab::HermitianMatrix<std::complex<double>, slab::lower> ml(3);
  for (unsigned i = 0; i < ml.n_rows(); ++i) {
    for (unsigned j = 0; j < i; ++j)
      ml(i, j) = std::complex<double>(3 * i + j, 3 * i + j);
    ml(i, i) = std::complex<double>(4 * i, 0);
  }
  std::cout << ml << std::endl;

  return 0;
}
