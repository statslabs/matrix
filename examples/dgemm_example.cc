//
// Created by Yi Pan (Institute of Cancer and Genomic Sciences) on 09/03/2018.
//

#include <cstdio>
#include <cstdlib>
#include "slab/matrix.h"

int min(int x, int y) { return (((x) < (y)) ? (x) : (y)); }

int main() {
  int m, n, k, i, j;

  printf(
      "\n This example computes real matrix C=alpha*A*B+beta*C using \n"
      " Intel(R) MKL function dgemm, where A, B, and  C are matrices and \n"
      " alpha and beta are double precision scalars\n\n");

  m = 2000, k = 200, n = 1000;
  printf(
      " Initializing data for matrix multiplication C=A*B for matrix \n"
      " A(%ix%i) and matrix B(%ix%i)\n\n",
      m, k, k, n);

  printf(
      " Allocating memory for matrices aligned on 64-byte boundary for better "
      "\n"
      " performance \n\n");
  slab::mat A(m, k);
  slab::mat B(k, n);
  slab::mat C(m, n);

  printf(" Intializing matrix data \n\n");
  for (i = 0; i < m; ++i) {
    for (j = 0; j < k; ++j) {
      A(i, j) = (double)(i * k + j + 1);
    }
  }

  for (i = 0; i < k; ++i) {
    for (j = 0; j < n; ++j) {
      B(i, j) = (double)(-(i * n + j + 1));
    }
  }

  printf(
      " Computing matrix product using Intel(R) MKL dgemm function via "
      "StatsLabs interface \n\n");

  C = slab::matmul(A, B);
  printf("\n Computations completed.\n\n");

  printf(" Top left corner of matrix A: \n");
  for (i = 0; i < min(m, 6); i++) {
    for (j = 0; j < min(k, 6); j++) {
      printf("%12.0f", A(i, j));
    }
    printf("\n");
  }

  printf("\n Top left corner of matrix B: \n");
  for (i = 0; i < min(k, 6); i++) {
    for (j = 0; j < min(n, 6); j++) {
      printf("%12.0f", B(i, j));
    }
    printf("\n");
  }

  printf("\n Top left corner of matrix C: \n");
  for (i = 0; i < min(m, 6); i++) {
    for (j = 0; j < min(n, 6); j++) {
      printf("%12.5G", C(i, j));
    }
    printf("\n");
  }

  printf(" Example completed. \n\n");
  return 0;
}
