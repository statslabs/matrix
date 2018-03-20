//
// Copyright 2018 The StatsLabs Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// -----------------------------------------------------------------------------
// blas.h
// -----------------------------------------------------------------------------
//

#ifndef SLAB_MATRIX_BLAS_H_
#define SLAB_MATRIX_BLAS_H_

#include "slab/matrix/matrix.h"
#include "slab/matrix/traits.h"

// Computes the parameters for a Givens rotation
//template<typename T>
//void blas_rotg(Matrix<T, 1> &a, Matrix<T, 1> &b, Matrix<T, 1> &c, Matrix<T, 1> &s) {
//  if (is_double<T>::value) {
//    cblas_drotg((double *) a.data(),
//                (double *) b.data(),
//                (double *) c.data(),
//                (double *) s.data());
//  } else if (is_float<T>::value) {
//    cblas_srotg((float *) a.data(), (float *) b.data(), (float *) c.data(), (float *) s.data());
//  } else if (is_complex_double<T>::value) {
//    cblas_zrotg((std::complex<double> *) a.data(), (const std::complex<double> *) b.data(),
//                (double *) c.data(), (std::complex<double> *) s.data());
//  } else if (is_complex_float<T>::value) {
//    cblas_crotg((std::complex<float> *) a.data(), (const std::complex<float> *) b.data(),
//                (float *) c.data(), (std::complex<float> *) s.data());
//  }
//}


// Swaps a vector with another vector
template<typename T>
void blas_swap(Matrix<T, 1> &x, Matrix<T, 1> &y)
{
  assert(x.size() == y.size());

  const int n = x.size();
  const int incx = 1;
  const int incy = 1;

  if (is_double<T>::value) {
    cblas_dswap(
        n,                   // n   : the number of elements in vectors x and y.
        (double *) x.data(), // x   : array, size at least (1 + (n-1)*abs(incx)).
        incx,                // incx: the increment for the elements of x.
        (double *) y.data(), // y   : array, size at least (1 + (n-1)*abs(incy))
        incy                 // incy: the increment for the elements of y.
    );
  } else if (is_float<T>::value) {
    cblas_sswap(n, (float *) x.data(), incx, (float *) y.data(), incy);
  } else if (is_complex_double<T>::value) {
    cblas_zswap(n, (std::complex<double> *) x.data(), incx, (std::complex<double> *) y.data(), incy);
  } else if (is_complex_float<T>::value) {
    cblas_cswap(n, (std::complex<float> *) x.data(), incx, (std::complex<float> *) y.data(), incy);
  }
}

// Computes the product of a vector by a scalar
template<typename T>
void blas_scal(const T a, Matrix<T, 1> &x)
{
  const int n = x.size();
  const int incx = 1;

  if (is_double<T>::value) {
    cblas_dscal(
        n,                   // n   : the number of elements in vector x
        (const double) a,    // a   : the scalar a.
        (double *) x.data(), // x   : array, size at least (1 + (n -1)*abs(incx))
        incx                 // incx: the increment for the elements of x.
    );
  } else if (is_float<T>::value) {
    cblas_sscal(n, (const float) a, (float *) x.data(), incx);
  } else if (is_complex_double<T>::value) {
    cblas_zdscal(n, (const double) a, (std::complex<double> *) x.data(), incx);
  } else if (is_complex_float<T>::value) {
    cblas_csscal(n, (const float) a, (std::complex<float> *) x.data(), incx);
  }
}

template<typename T>
void blas_scal(const std::complex<T> &a, Matrix<std::complex<T>, 1> &x)
{
  const int n = x.size();
  const int incx = 1;

  if (is_double<T>::value) {
    cblas_zscal(
        n,
        (const std::complex<double> *) &a,
        (std::complex<double> *) x.data(),
        incx
    );
  } else if (is_float<T>::value) {
    cblas_cscal(
        n,
        (const std::complex<float> *) &a,
        (std::complex<float> *) x.data(),
        incx);
  }
}


// Copies vector to another vector
template<typename T>
void blas_copy(const Matrix<T, 1> &x, Matrix<T, 1> &y) {
  y.clear();
  y = Matrix<T, 1>(x.size());

  const int incx = 1;
  const int incy = 1;

  if (is_double<T>::value)
    cblas_dcopy(
        x.size(),
        (const double *) x.data(),
        incx,
        (double *) y.data(),
        incy
    );
  else if (is_float<T>::value) {
    cblas_scopy(x.size(), (const float *) x.data(), incx, (float *) y.data(), incy);
  } else if (is_complex_double<T>::value) {
    cblas_zcopy(x.size(), (const std::complex<double> *) x.data(), incx, (std::complex<double> *) y.data(), incy);
  } else if (is_complex_float<T>::value) {
    cblas_ccopy(x.size(), (const std::complex<float> *) x.data(), incx, (std::complex<float> *) y.data(), incy);
  }
}

// Computes a vector-scalar product and adds the result to a vector
template<typename T>
void blas_axpy(const T &a, const Matrix<T, 1> &x, Matrix<T, 1> &y) {
  assert(x.size() == y.size());

  const int n = x.size();
  const int incx = 1;
  const int incy = 1;

  if (is_double<T>::value) {
    cblas_daxpy(
        n,                         // n   : the number of elements in vectors x and y.
        a,                         // a   : the scalar a.
        (const double *) x.data(), // x   : the vector x.
        incx,                      // incx: the increment for indexing vector x.
        (double *) y.data(),       // y   : the vector y.
        incy                       // incy: the increment for indexing vector y.
    );
  } else if (is_float<T>::value) {
    cblas_saxpy(n, a, (const float *) x.data(), incx, (float *) y.data(), incy);
  } else if (is_complex_double<T>::value) {
    cblas_zaxpy(n, &a, (const std::complex<double> *) x.data(), incx, (std::complex<double> *) y.data(), incy);
  } else if (is_complex_float<T>::value) {
    cblas_caxpy(n, &a, (const std::complex<float> *) x.data(), incx, (std::complex<float> *) y.data(), incy);
  }
}

// Computes the Euclidean norm of a vector
template<typename T>
double blas_nrm2(const Matrix<T, 1> &x) {
  double res = 0.0;

  const int n = x.size();
  const int incx = 1;

  if (is_double<T>::value) {
    res = cblas_dnrm2(
        n,
        (const double *) x.data(),
        incx
    );
  } else if (is_float<T>::value) {
    res = cblas_snrm2(n, (const float *) x.data(), incx);
  } else if (is_complex_double<T>::value) {
    res = cblas_dznrm2(n, (const std::complex<double> *) x.data(), incx);
  } else if (is_complex_float<T>::value) {
    res = cblas_scnrm2(n, (const std::complex<float> *) x.data(), incx);
  }

  return res;
}

// Computes the sum of magnitudes of the vector elements
template<typename T>
T blas_asum(const Matrix<T, 1> &x) {
  const int n = x.size();
  const int incx = 1;

  T res;
  if (is_double<T>::value) {
    res = cblas_dasum(
        n,                         // n   : the number of elements in vector x.
        (const double *) x.data(), // x   : the vector x.
        incx                       // incx: the increment for indexing vector x.
    );
  } else if (is_float<T>::value) {
    res = cblas_sasum(n, (const float *) x.data(), incx);
  } else if (is_complex_double<T>::value) {
    res = cblas_dzasum(n, (const std::complex<double> *) x.data(), incx);
  } else if (is_complex_float<T>::value) {
    res = cblas_scasum(n, (const std::complex<float> *) x.data(), incx);
  }

  return res;
}

// Finds the index of the element with maximum absolute value
template<typename T>
std::size_t blas_iamax(const Matrix<T, 1> &x) {
  std::size_t res = 0;
  std::size_t incx = 1;

  if (is_double<T>::value) {
    res = cblas_idamax(
        x.size(),
        (const double *) x.data(),
        incx
    );
  } else if (is_float<T>::value) {
    res = cblas_isamax(x.size(), (const float *) x.data(), incx);
  } else if (is_complex_double<T>::value) {
    res = cblas_izamax(x.size(), (const std::complex<double> *) x.data(), incx);
  } else if (is_complex_float<T>::value) {
    res = cblas_icamax(x.size(), (const std::complex<float> *) x.data(), incx);
  }

  return res;
}

#endif // SLAB_MATRIX_BLAS_H_
