#include <exception>
#include <iostream>
using namespace std;

#include "slab/matrix.h"
using namespace slab;

using Mat2d = Matrix<double, 2>;
using Vec = Matrix<double, 1>;

struct Elim_failure : public exception {
  int row;
  Elim_failure(int r) : row(r) {}
  const char *what() const noexcept { return "elimination failure"; }
};

struct Back_subst_failure : public exception {
  int row;
  Back_subst_failure(int r) : row(r) {}
  const char *what() const noexcept { return "back substitution failure"; }
};

Vec scale_and_add(const MatrixRef<double, 1> &v1, double s,
                  const MatrixRef<double, 1> &v2) {
  Vec res = v2;
  blas_axpy(s, v1, res);

  return res;
}

template <typename T, typename T1>
double dot_product(const MatrixRef<T, 1> &a, const MatrixRef<T1, 1> &b) {
  return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

void classical_elimination(Mat2d &A, Vec &b) {
  const std::size_t n = A.n_rows();

  // traverse from 1st column to the next-to-last, filling zeros into all
  // elements under the diagonal:
  for (std::size_t j = 0; j != n - 1; ++j) {
    const double pivot = A(j, j);
    if (pivot == 0) throw Elim_failure(j);
    // fill zeros into each element under the diagonal of the ith row:
    for (std::size_t i = j + 1; i != n; ++i) {
      const double mult = A(i, j) / pivot;
      A[i](slice(j)) = scale_and_add(A[j](slice(j)), -mult, A[i](slice(j)));
      b(i) -= mult * b(j);  // make the corresponding change to b
    }
  }
}

Vec back_substitution(const Mat2d &A, const Vec &b) {
  const std::size_t n = A.n_rows();
  Vec x(n);

  for (int i = n - 1; i >= 0; --i) {
    double s = b(i) - dot_product(A[i](slice(i + 1)), x(slice(i + 1)));
    if (double m = A(i, i))
      x(i) = s / m;
    else
      throw Back_subst_failure(i);
  }
  return x;
}

Vec classical_gaussian_elimination(Mat2d A, Vec b) {
  classical_elimination(A, b);
  return back_substitution(A, b);
}

int main() {
  Mat2d A = {{2, 1, -1}, {-3, -1, 2}, {-2, 1, 2}};

  Vec b = {8, -11, -3};

  Vec res = classical_gaussian_elimination(A, b);

  cout << "res = " << res << endl;

  return 0;
}