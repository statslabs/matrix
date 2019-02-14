#ifndef STATSLABS_MATRIX_TEST_MATRIX_FNS_H_
#define STATSLABS_MATRIX_TEST_MATRIX_FNS_H_

#include "slab/matrix.h"

namespace slab {

TEST(MatrixFunctions, Prod) {
  vec v = {1, 2, 3};
  EXPECT_EQ(6, prod(v));

  mat m = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9}};
  vec m_prod_ans = {28, 80, 162};
  EXPECT_EQ(m_prod_ans, prod(m));
}

} // namespace slab

#endif  // STATSLABS_MATRIX_TEST_MATRIX_FNS_H_
