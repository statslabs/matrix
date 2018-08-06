#ifndef SLAB_MATRIX_TYPEDEF_H_
#define SLAB_MATRIX_TYPEDEF_H_

// General Matrix -- Vector / Matrix / Cube

using vec      = Matrix<double, 1>;
using mat      = Matrix<double, 2>;
using cube     = Matrix<double, 3>;

using dvec     = Matrix<double, 1>;
using dmat     = Matrix<double, 2>;
using dcube    = Matrix<double, 3>;

using fvec     = Matrix<float, 1>;
using fmat     = Matrix<float, 2>;
using fcube    = Matrix<float, 3>;

using cx_vec   = Matrix<std::complex<double>, 1>;
using cx_mat   = Matrix<std::complex<double>, 2>;
using cx_cube  = Matrix<std::complex<double>, 3>;

using cx_dvec  = Matrix<std::complex<double>, 1>;
using cx_dmat  = Matrix<std::complex<double>, 2>;
using cx_dcube = Matrix<std::complex<double>, 3>;

using cx_fvec  = Matrix<std::complex<float>, 1>;
using cx_fmat  = Matrix<std::complex<float>, 2>;
using cx_fcube = Matrix<std::complex<float>, 3>;

using uvec     = Matrix<unsigned, 1>;
using umat     = Matrix<unsigned, 2>;
using ucube    = Matrix<unsigned, 3>;

using ivec     = Matrix<int, 1>;
using imat     = Matrix<int, 2>;
using icube    = Matrix<int, 3>;

// Packed Matrix -- Symmetric Matrix / Triangular Matrix / Hermitian Matrix

using symm_mat     = SymmetricMatrix<double, upper>;
using utri_mat     = TriangularMatrix<double, upper>;
using ltri_mat     = TriangularMatrix<double, lower>;
using herm_mat     = HermitianMatrix<double, upper>;

using symm_dmat    = SymmetricMatrix<double, upper>;
using utri_dmat    = TriangularMatrix<double, upper>;
using ltri_dmat    = TriangularMatrix<double, lower>;
using herm_dmat    = HermitianMatrix<double, upper>;

using symm_fmat    = SymmetricMatrix<float, upper>;
using utri_fmat    = TriangularMatrix<float, upper>;
using ltri_fmat    = TriangularMatrix<float, lower>;
using herm_fmat    = HermitianMatrix<float, upper>;

using symm_cx_mat  = SymmetricMatrix<std::complex<double>, upper>;
using utri_cx_mat  = TriangularMatrix<std::complex<double>, upper>;
using ltri_cx_mat  = TriangularMatrix<std::complex<double>, lower>;
using herm_cx_mat  = HermitianMatrix<std::complex<double>, upper>;

using symm_cx_dmat = SymmetricMatrix<std::complex<double>, upper>;
using utri_cx_dmat = TriangularMatrix<std::complex<double>, upper>;
using ltri_cx_dmat = TriangularMatrix<std::complex<double>, lower>;
using herm_cx_dmat = HermitianMatrix<std::complex<double>, upper>;

using symm_cx_fmat = SymmetricMatrix<std::complex<float>, upper>;
using utri_cx_fmat = TriangularMatrix<std::complex<float>, upper>;
using ltri_cx_fmat = TriangularMatrix<std::complex<float>, lower>;
using herm_cx_fmat = HermitianMatrix<std::complex<float>, upper>;

#endif // SLAB_MATRIX_TYPEDEF_H_
