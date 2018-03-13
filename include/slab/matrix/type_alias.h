#ifndef SLAB_MATRIX_TYPEDEF_H_
#define SLAB_MATRIX_TYPEDEF_H_

#include "slab/matrix/matrix.h"

using vec   = Matrix<double, 1>;
using mat   = Matrix<double, 2>;
using cube  = Matrix<double, 3>;

using fvec  = Matrix<float, 1>;
using fmat  = Matrix<float, 2>;
using fcube = Matrix<float, 3>;

using dvec  = Matrix<double, 1>;
using dmat  = Matrix<double, 2>;
using dcube = Matrix<double, 3>;

#endif // SLAB_MATRIX_TYPEDEF_H_
