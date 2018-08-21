/*
   LAPACKE_dgels Example.
   ======================

   Program computes the least squares solution to the overdetermined linear
   system A*X = B with full rank matrix A using QR factorization,
   where A is the coefficient matrix:

     1.44  -7.84  -4.39   4.53
    -9.96  -0.28  -3.24   3.83
    -7.55   3.24   6.27  -6.64
     8.34   8.09   5.28   2.06
     7.08   2.52   0.74  -2.47
    -5.45  -5.70  -1.19   4.70

   and B is the right-hand side matrix:

     8.58   9.35
     8.26  -4.43
     8.48  -0.70
    -5.28  -0.26
     5.72  -7.36
     8.93  -2.52

   Description.
   ============

   The routine solves overdetermined or underdetermined real linear systems
   involving an m-by-n matrix A, or its transpose, using a QR or LQ
   factorization of A. It is assumed that A has full rank.

   Several right hand side vectors b and solution vectors x can be handled
   in a single call; they are stored as the columns of the m-by-nrhs right
   hand side matrix B and the n-by-nrhs solution matrix X.

   Example Program Results.
   ========================

 LAPACKE_dgels (row-major, high-level) Example Program Results

 Solution
  -0.45   0.25
  -0.85  -0.90
   0.71   0.63
   0.13   0.14

 Residual sum of squares for the solution
 195.36 107.06

 Details of QR factorization
 -17.54  -4.76  -1.96   0.42
  -0.52  12.40   7.88  -5.84
  -0.40  -0.14  -5.75   4.11
   0.44  -0.66  -0.20  -7.78
   0.37  -0.26  -0.17  -0.15
  -0.29   0.46   0.41   0.24
*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "slab/matrix.h"

/* Main program */
int main() {
  /* Locals */
  int info;
  /* Local arrays */
  slab::mat a = {{1.44, -7.84, -4.39, 4.53}, {-9.96, -0.28, -3.24, 3.83},
                 {-7.55, 3.24, 6.27, -6.64}, {8.34, 8.09, 5.28, 2.06},
                 {7.08, 2.52, 0.74, -2.47},  {-5.45, -5.70, -1.19, 4.70}};
  slab::mat b = {{8.58, 9.35},   {8.26, -4.43}, {8.48, -0.70},
                 {-5.28, -0.26}, {5.72, -7.36}, {8.93, -2.52}};
  /* Executable statements */
  printf("lapack_dgels Example Program Results\n");
  /* Solve the equations AX = B */
  info = slab::lapack_gels('N', a, b);
  /* Check for the full rank */
  if (info > 0) {
    printf("The diagonal element %i of the triangular factor ", info);
    printf("of A is zero, so that A does not have full rank;\n");
    printf("the least squares solution could not be computed.\n");
    exit(1);
  }
  /* Print least squares solution */
  b.rows(0, a.n_cols() - 1).print("Least squares solution");
  /* Print residual sum of squares for the solution */
  printf("\n %s\n", "Residual sum of squares for the solution");
  for (std::size_t j = 0; j < b.n_cols(); ++j) {
    slab::vec tmp =
        b(slab::slice{a.n_cols(), a.n_rows() - a.n_cols()}, slab::slice{j, 1});
    double norm = slab::dot(tmp, tmp);
    printf(" %6.2f", norm);
  }
  printf("\n");
  /* Print details of QR factorization */
  a.print("Details of QR factorization");
  exit(0);
} /* End of lapack_gels Example */
