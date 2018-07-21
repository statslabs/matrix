//
// Copyright 2018 The Statslabs Authors.
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

/// @file blas_interface.h
/// @brief BLAS interface

#ifndef SLAB_MATRIX_BLAS_INTERFACE_H_
#define SLAB_MATRIX_BLAS_INTERFACE_H_

#include "slab/matrix/matrix.h"
#include "slab/matrix/traits.h"

enum class blas_trans {
  no_trans = CblasNoTrans,
  trans = CblasTrans,
  conj_trans = CblasConjTrans
};

/// @addtogroup blas_interface BLAS INTERFACE
/// @{

/// @addtogroup blas_level1 BLAS Level 1
/// @{

/// @brief Computes the sum of magnitudes of the vector elements
/// @param x a vector.
template<typename T>
T blas_asum(const Matrix<T, 1> &x) {
  const int n = x.size();
  const int incx = x.descriptor().strides[0];

  T res;
  if (is_double<T>::value) {
    res = cblas_dasum(
        n,
        (const double *) x.data(),
        incx
    );
  } else if (is_float<T>::value) {
    res = cblas_sasum(
        n,
        (const float *) x.data(),
        incx
    );
  } else if (is_complex_double<T>::value) {
    res = cblas_dzasum(
        n,
        (const std::complex<double> *) x.data(),
        incx
    );
  } else if (is_complex_float<T>::value) {
    res = cblas_scasum(
        n,
        (const std::complex<float> *) x.data(),
        incx
    );
  }

  return res;
}

/// @brief Computes a vector-scalar product and adds the result to a vector
template<typename T>
void blas_axpy(const T &a, const MatrixBase<T, 1> &x, MatrixBase<T, 1> &y) {
  assert(x.size() == y.size());

  const int n = x.size();
  const int incx = x.descriptor().strides[0];
  const int incy = y.descriptor().strides[0];

  if (is_double<T>::value) {
    cblas_daxpy(
        n,
        a,
        (const double *) (x.data() + x.descriptor().start),
        incx,
        (double *) (y.data() + y.descriptor().start),
        incy
    );
  } else if (is_float<T>::value) {
    cblas_saxpy(
        n,
        a,
        (const float *) (x.data() + x.descriptor().start),
        incx,
        (float *) (y.data() + y.descriptor().start),
        incy
    );
  } else if (is_complex_double<T>::value) {
    cblas_zaxpy(
        n,
        &a,
        (const std::complex<double> *) (x.data() + x.descriptor().start),
        incx,
        (std::complex<double> *) (y.data() + y.descriptor().start),
        incy
    );
  } else if (is_complex_float<T>::value) {
    cblas_caxpy(
        n,
        &a,
        (const std::complex<float> *) (x.data() + x.descriptor().start),
        incx,
        (std::complex<float> *) (y.data() + y.descriptor().start),
        incy
    );
  }
}

/// @brief Copies vector to another vector
template<typename T>
void blas_copy(const Matrix<T, 1> &x, Matrix<T, 1> &y) {
  y.clear();
  y = Matrix<T, 1>(x.size());

  const int incx = x.descriptor().strides[0];
  const int incy = y.descriptor().strides[0];

  if (is_double<T>::value)
    cblas_dcopy(
        x.size(),
        (const double *) (x.data() + x.descriptor().start),
        incx,
        (double *) (y.data() + y.descriptor().start),
        incy
    );
  else if (is_float<T>::value) {
    cblas_scopy(
        x.size(),
        (const float *) (x.data() + x.descriptor().start),
        incx,
        (float *) (y.data() + y.descriptor().start),
        incy
    );
  } else if (is_complex_double<T>::value) {
    cblas_zcopy(
        x.size(),
        (const std::complex<double> *) (x.data() + x.descriptor().start),
        incx,
        (std::complex<double> *) (y.data() + y.descriptor().start),
        incy
    );
  } else if (is_complex_float<T>::value) {
    cblas_ccopy(
        x.size(),
        (const std::complex<float> *) (x.data() + x.descriptor().start),
        incx,
        (std::complex<float> *) (y.data() + y.descriptor().start),
        incy
    );
  }
}

/// @brief Computes a vector-vector dot product
template<typename T>
T blas_dot(const Matrix<T, 1> &x, const Matrix<T, 1> &y) {
  assert(x.size() == y.size());

  const int n = x.size();
  const int incx = x.descriptor().strides[0];
  const int incy = y.descriptor().strides[0];

  T res = 0.0;
  if (is_double<T>::value) {
    res = cblas_ddot(
        n,
        (const double *) (x.data() + x.descriptor().start),
        incx,
        (const double *) (y.data() + y.descriptor().start),
        incy
    );
  } else if (is_float<T>::value) {
    res = cblas_sdot(
        n,
        (const float *) (x.data() + x.descriptor().start),
        incx,
        (const float *) (y.data() + y.descriptor().start),
        incy
    );
  }

  return res;
}

/// @brief Computes a vector-vector dot product with double precision
float blas_sdsdot(const float sb, const Matrix<float, 1> &sx,
                  const Matrix<float, 1> &sy) {
  assert(sx.size() == sy.size());

  const int n = sx.size();
  const int incx = sx.descriptor().strides[0];
  const int incy = sy.descriptor().strides[0];

  return cblas_sdsdot(n, sb, sx.data() + sx.descriptor().start, incx,
                      sy.data() + sy.descriptor().start, incy);
}

/// @brief Computes a vector-vector dot product with double precision
double blas_dsdot(const Matrix<float, 1> &sx, const Matrix<float, 1> &sy) {
  assert(sx.size() == sy.size());

  const int n = sx.size();
  const int incx = sx.descriptor().strides[0];
  const int incy = sy.descriptor().strides[0];

  return cblas_dsdot(n, sx.data() + sx.descriptor().start, incx,
                     sy.data() + sy.descriptor().start, incy);
}

/// @brief Computes the Euclidean norm of a vector
template<typename T>
double blas_nrm2(const Matrix<T, 1> &x) {
  double res = 0.0;

  const int n = x.size();
  const int incx = x.descriptor().strides[0];

  if (is_double<T>::value) {
    res = cblas_dnrm2(
        n,
        (const double *) (x.data() + x.descriptor().start),
        incx
    );
  } else if (is_float<T>::value) {
    res = cblas_snrm2(
        n,
        (const float *) (x.data() + x.descriptor().start),
        incx
    );
  } else if (is_complex_double<T>::value) {
    res = cblas_dznrm2(
        n,
        (const std::complex<double> *) (x.data() + x.descriptor().start),
        incx
    );
  } else if (is_complex_float<T>::value) {
    res = cblas_scnrm2(
        n,
        (const std::complex<float> *) (x.data() + x.descriptor().start),
        incx
    );
  }

  return res;
}

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

/// @brief Computes the product of a vector by a scalar
template<typename T>
void blas_scal(const T a, Matrix<T, 1> &x) {
  const int n = x.size();
  const int incx = x.descriptor().strides[0];

  if (is_double<T>::value) {
    cblas_dscal(
        n,
        (const double) a,
        (double *) (x.data() + x.descriptor().start),
        incx
    );
  } else if (is_float<T>::value) {
    cblas_sscal(
        n,
        (const float) a,
        (float *) (x.data() + x.descriptor().start),
        incx
    );
  } else if (is_complex_double<T>::value) {
    cblas_zdscal(
        n,
        (const double) a,
        (std::complex<double> *) (x.data() + x.descriptor().start),
        incx
    );
  } else if (is_complex_float<T>::value) {
    cblas_csscal(
        n,
        (const float) a,
        (std::complex<float> *) (x.data() + x.descriptor().start),
        incx
    );
  }
}

template<typename T>
void blas_scal(const std::complex<T> &a, Matrix<std::complex<T>, 1> &x) {
  const int n = x.size();
  const int incx = x.descriptor().strides[0];

  if (is_double<T>::value) {
    cblas_zscal(
        n,
        (const std::complex<double> *) &a,
        (std::complex<double> *) (x.data() + x.descriptor().start),
        incx
    );
  } else if (is_float<T>::value) {
    cblas_cscal(
        n,
        (const std::complex<float> *) &a,
        (std::complex<float> *) (x.data() + x.descriptor().start),
        incx);
  }
}

/// @brief Swaps a vector with another vector.
///
/// @param x a vector.
/// @param y another vector.
template<typename T>
void blas_swap(Matrix<T, 1> &x, Matrix<T, 1> &y) {
  assert(x.size() == y.size());

  const int n = x.size();
  const int incx = x.descriptor().strides[0];
  const int incy = y.descriptor().strides[0];

  if (is_double<T>::value) {
    cblas_dswap(
        n,
        (double *) (x.data() + x.descriptor().start),
        incx,
        (double *) (y.data() + y.descriptor().start),
        incy
    );
  } else if (is_float<T>::value) {
    cblas_sswap(
        n,
        (float *) (x.data() + x.descriptor().start),
        incx,
		(float *) (y.data() + y.descriptor().start),
        incy
    );
  } else if (is_complex_double<T>::value) {
    cblas_zswap(
        n,
        (std::complex<double> *) (x.data() + x.descriptor().start),
		incx,
        (std::complex<double> *) (y.data() + y.descriptor().start),
        incy
    );
  } else if (is_complex_float<T>::value) {
    cblas_cswap(
        n,
        (std::complex<float> *) (x.data() + x.descriptor().start),
        incx,
        (std::complex<float> *) (y.data() + y.descriptor().start),
        incy
    );
  }
}

/// @brief Finds the index of the element with maximum absolute value
template<typename T>
std::size_t blas_iamax(const Matrix<T, 1> &x) {
  std::size_t res = 0;
  std::size_t incx = x.descriptor().strides[0];

  if (is_double<T>::value) {
    res = cblas_idamax(
        x.size(),
        (const double *) (x.data() + x.descriptor().start),
        incx
    );
  } else if (is_float<T>::value) {
    res = cblas_isamax(
        x.size(),
        (const float *) (x.data() + x.descriptor().start),
        incx
    );
  } else if (is_complex_double<T>::value) {
    res = cblas_izamax(
        x.size(),
        (const std::complex<double> *) (x.data() + x.descriptor().start),
        incx
    );
  } else if (is_complex_float<T>::value) {
    res = cblas_icamax(
        x.size(),
        (const std::complex<float> *) (x.data() + x.descriptor().start),
        incx
    );
  }

  return res;
}
/// @}

/// @addtogroup blas_level2 BLAS Level 2
/// @{

void blas_gemv()
{

}

/// @}

/// @addtogroup blas_level3 BLAS Level 3
/// @{

void blas_gemm()
{

}

/// @}
/// @} BLAS INTERFACE

#endif // SLAB_MATRIX_BLAS_INTERFACE_H_
