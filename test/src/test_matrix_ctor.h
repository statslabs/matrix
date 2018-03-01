#ifndef MATRIX_TEST_MATRIX_CTOR_H_
#define MATRIX_TEST_MATRIX_CTOR_H_

#include <array>
#include <gtest/gtest.h>
#include "slab/matrix.h"

namespace slab {

TEST(MatrixCtorTest, ConstructFromMatrixRef) {
  Matrix<double, 1> m1{1, 2, 3};
  auto mr1 = m1(slice(1));                    // MatrixRef<double, 1>
  Matrix<double, 1> m1sub(mr1);

  EXPECT_EQ(2, m1sub(0));
  EXPECT_EQ(3, m1sub(1));

  Matrix<double, 2> m2{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  auto mr2 = m2(slice(1), slice(1));          // MatrixRef<double, 2>
  Matrix<double, 2> m2sub(mr2);

  EXPECT_EQ(5, m2sub(0, 0));
  EXPECT_EQ(6, m2sub(0, 1));
  EXPECT_EQ(8, m2sub(1, 0));
  EXPECT_EQ(9, m2sub(1, 1));

  Matrix<double, 3>
      m3{{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
         {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
         {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}};
  auto mr3 = m3(slice(1), slice(1), slice(1)); // MatrixRef<double, 3>
  Matrix<double, 3> m3sub(mr3);

  EXPECT_EQ(5, m3sub(0, 0, 0));
  EXPECT_EQ(6, m3sub(0, 0, 1));
  EXPECT_EQ(8, m3sub(0, 1, 0));
  EXPECT_EQ(9, m3sub(0, 1, 1));
  EXPECT_EQ(5, m3sub(1, 0, 0));
  EXPECT_EQ(6, m3sub(1, 0, 1));
  EXPECT_EQ(8, m3sub(1, 1, 0));
  EXPECT_EQ(9, m3sub(1, 1, 1));
}

TEST(MatrixCtorTest, ConstructFromExtent) {

  Matrix<double, 0> m0;
  auto ms0 = m0.descriptor(); // MatrixSlice<0>

  EXPECT_EQ(m0(), 0);         // the only element
  EXPECT_EQ(ms0.size, 0);
  EXPECT_EQ(ms0.start, 0);

  std::array<int, 3> arr = {3, 4, 5};

  Matrix<double, 1> m1(3);
  auto ms1 = m1.descriptor(); // MatrixSlice<1>

  EXPECT_EQ(m1(0), 0);        // first element
  EXPECT_EQ(m1(2), 0);        // last element
  EXPECT_EQ(ms1.size, 3);
  EXPECT_EQ(ms1.start, 0);
  for (auto idx = 0; idx < m1.order; ++idx) {
    EXPECT_EQ(ms1.extents[idx], arr[idx]) << "arrays differ at index " << idx;
  }

  Matrix<double, 2> m2(3, 4);
  auto ms2 = m2.descriptor(); // MatrixSlice<2>

  EXPECT_EQ(m2(0, 0), 0);     // first element
  EXPECT_EQ(m2(2, 3), 0);     // last element
  EXPECT_EQ(ms2.size, 12);
  EXPECT_EQ(ms2.start, 0);
  for (auto idx = 0; idx < m2.order; ++idx) {
    EXPECT_EQ(ms2.extents[idx], arr[idx]) << "arrays differ at index " << idx;
  }

  Matrix<double, 3> m3(3, 4, 5);
  auto ms3 = m3.descriptor(); // MatrixSlice<3>

  EXPECT_EQ(m3(0, 0, 0), 0);  // first element
  EXPECT_EQ(m3(2, 3, 4), 0);  // last element
  EXPECT_EQ(ms3.size, 60);
  EXPECT_EQ(ms3.start, 0);
  for (auto idx = 0; idx < m3.order; ++idx) {
    EXPECT_EQ(ms3.extents[idx], arr[idx]) << "arrays differ at index " << idx;
  }
}

TEST(MatrixCtorTest, ConstructFromMatrixInitializer) {
  Matrix<double, 0> m0{1};

  EXPECT_EQ(m0(), 1);

  Matrix<double, 1> m1{1, 2};

  EXPECT_EQ(m1(0), 1);
  EXPECT_EQ(m1(1), 2);

  Matrix<double, 2> m2{{1, 2}, {3, 4}};

  EXPECT_EQ(m2(0, 0), 1);
  EXPECT_EQ(m2(0, 1), 2);
  EXPECT_EQ(m2(1, 0), 3);
  EXPECT_EQ(m2(1, 1), 4);

  Matrix<double, 3> m3{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}};

  EXPECT_EQ(m3(0, 0, 0), 1);
  EXPECT_EQ(m3(0, 0, 1), 2);
  EXPECT_EQ(m3(0, 1, 0), 3);
  EXPECT_EQ(m3(0, 1, 1), 4);
  EXPECT_EQ(m3(1, 0, 0), 5);
  EXPECT_EQ(m3(1, 0, 1), 6);
  EXPECT_EQ(m3(1, 1, 0), 7);
  EXPECT_EQ(m3(1, 1, 1), 8);
}

} // namespace slab

#endif // MATRIX_TEST_MATRIX_CTOR_H_
