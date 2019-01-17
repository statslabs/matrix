#ifndef MATRIX_TEST_CONSTRUCT_AND_ASSIGNMENT_H
#define MATRIX_TEST_CONSTRUCT_AND_ASSIGNMENT_H

#include <array>

#include <gtest/gtest.h>
#include "slab/matrix.h"

namespace slab {

TEST(MatrixConstructionTest, ConstructFromExtent) {
  Matrix<double, 0> m0;
  auto ms0 = m0.descriptor();  // MatrixSlice<0>

  EXPECT_EQ(0, m0());  // the only element
  EXPECT_EQ(0, ms0.size);
  EXPECT_EQ(0, ms0.start);

  std::array<int, 3> arr = {{3, 4, 5}};

  Matrix<double, 1> m1(3);
  auto ms1 = m1.descriptor();  // MatrixSlice<1>

  EXPECT_EQ(0, m1(0));  // first element
  EXPECT_EQ(0, m1(2));  // last element
  EXPECT_EQ(3, ms1.size);
  EXPECT_EQ(0, ms1.start);
  for (std::size_t idx = 0; idx < m1.order(); ++idx) {
    EXPECT_EQ(arr[idx], ms1.extents[idx]) << "arrays differ at index " << idx;
  }

  Matrix<double, 2> m2(3, 4);
  auto ms2 = m2.descriptor();  // MatrixSlice<2>

  EXPECT_EQ(0, m2(0, 0));  // first element
  EXPECT_EQ(0, m2(2, 3));  // last element
  EXPECT_EQ(12, ms2.size);
  EXPECT_EQ(0, ms2.start);
  for (std::size_t idx = 0; idx < m2.order(); ++idx) {
    EXPECT_EQ(arr[idx], ms2.extents[idx]) << "arrays differ at index " << idx;
  }

  Matrix<double, 3> m3(3, 4, 5);
  auto ms3 = m3.descriptor();  // MatrixSlice<3>

  EXPECT_EQ(0, m3(0, 0, 0));  // first element
  EXPECT_EQ(0, m3(2, 3, 4));  // last element
  EXPECT_EQ(60, ms3.size);
  EXPECT_EQ(0, ms3.start);
  for (std::size_t idx = 0; idx < m3.order(); ++idx) {
    EXPECT_EQ(arr[idx], ms3.extents[idx]) << "arrays differ at index " << idx;
  }
}

TEST(MatrixConstructionTest, ConstructFromMatrixInitializer) {
  Matrix<double, 0> m0{1};

  EXPECT_EQ(1, m0());

  Matrix<double, 1> m1{1, 2};

  EXPECT_EQ(1, m1(0));
  EXPECT_EQ(2, m1(1));

  Matrix<double, 2> m2{{1, 2}, {3, 4}};

  EXPECT_EQ(1, m2(0, 0));
  EXPECT_EQ(2, m2(0, 1));
  EXPECT_EQ(3, m2(1, 0));
  EXPECT_EQ(4, m2(1, 1));

  Matrix<double, 3> m3{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}};

  EXPECT_EQ(1, m3(0, 0, 0));
  EXPECT_EQ(2, m3(0, 0, 1));
  EXPECT_EQ(3, m3(0, 1, 0));
  EXPECT_EQ(4, m3(0, 1, 1));
  EXPECT_EQ(5, m3(1, 0, 0));
  EXPECT_EQ(6, m3(1, 0, 1));
  EXPECT_EQ(7, m3(1, 1, 0));
  EXPECT_EQ(8, m3(1, 1, 1));
}

TEST(MatrixConstructionTest, ConstructFromMatrixRef) {
  Matrix<double, 1> m1{1, 2, 3};
  auto mr1 = m1(slice(1));  // MatrixRef<double, 1>
  Matrix<double, 1> m1sub(mr1);

  EXPECT_EQ(2, m1sub(0));
  EXPECT_EQ(3, m1sub(1));

  Matrix<double, 2> m2{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  auto mr2 = m2(slice(1), slice(1));  // MatrixRef<double, 2>
  Matrix<double, 2> m2sub(mr2);

  EXPECT_EQ(5, m2sub(0, 0));
  EXPECT_EQ(6, m2sub(0, 1));
  EXPECT_EQ(8, m2sub(1, 0));
  EXPECT_EQ(9, m2sub(1, 1));

  Matrix<double, 3> m3{{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
                       {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
                       {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}};
  auto mr3 = m3(slice(1), slice(1), slice(1));  // MatrixRef<double, 3>
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

TEST(MatrixAssignmentTest, AssignFromMatrixInitializer) {
  Matrix<double, 0> m0 = {1};

  EXPECT_EQ(1, m0());

  Matrix<double, 1> m1 = {1, 2};

  EXPECT_EQ(1, m1(0));
  EXPECT_EQ(2, m1(1));

  Matrix<double, 2> m2 = {{1, 2}, {3, 4}};

  EXPECT_EQ(1, m2(0, 0));
  EXPECT_EQ(2, m2(0, 1));
  EXPECT_EQ(3, m2(1, 0));
  EXPECT_EQ(4, m2(1, 1));

  Matrix<double, 3> m3 = {{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}};

  EXPECT_EQ(1, m3(0, 0, 0));
  EXPECT_EQ(2, m3(0, 0, 1));
  EXPECT_EQ(3, m3(0, 1, 0));
  EXPECT_EQ(4, m3(0, 1, 1));
  EXPECT_EQ(5, m3(1, 0, 0));
  EXPECT_EQ(6, m3(1, 0, 1));
  EXPECT_EQ(7, m3(1, 1, 0));
  EXPECT_EQ(8, m3(1, 1, 1));
}

TEST(MatrixAssignmentTest, AssignFromMatrixRef) {
  Matrix<double, 1> m1{1, 2, 3};
  auto mr1 = m1(slice(1));  // MatrixRef<double, 1>
  Matrix<double, 1> m1sub = mr1;

  EXPECT_EQ(2, m1sub(0));
  EXPECT_EQ(3, m1sub(1));

  Matrix<double, 2> m2{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  auto mr2 = m2(slice(1), slice(1));  // MatrixRef<double, 2>
  Matrix<double, 2> m2sub = mr2;

  EXPECT_EQ(5, m2sub(0, 0));
  EXPECT_EQ(6, m2sub(0, 1));
  EXPECT_EQ(8, m2sub(1, 0));
  EXPECT_EQ(9, m2sub(1, 1));

  Matrix<double, 3> m3{{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
                       {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
                       {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}};
  auto mr3 = m3(slice(1), slice(1), slice(1));  // MatrixRef<double, 3>
  Matrix<double, 3> m3sub = mr3;

  EXPECT_EQ(5, m3sub(0, 0, 0));
  EXPECT_EQ(6, m3sub(0, 0, 1));
  EXPECT_EQ(8, m3sub(0, 1, 0));
  EXPECT_EQ(9, m3sub(0, 1, 1));
  EXPECT_EQ(5, m3sub(1, 0, 0));
  EXPECT_EQ(6, m3sub(1, 0, 1));
  EXPECT_EQ(8, m3sub(1, 1, 0));
  EXPECT_EQ(9, m3sub(1, 1, 1));
}

}  // namespace slab

#endif  // MATRIX_TEST_CONSTRUCT_AND_ASSIGNMENT_H
