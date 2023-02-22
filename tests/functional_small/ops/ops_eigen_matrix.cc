
#include <gtest/gtest.h>
#include "pressio/ops.hpp"

TEST(ops_eigen_matrix, dense_matrix_clone)
{
  using T = Eigen::MatrixXd;
  T a(6,8);

  int c=0;
  for (int i=0; i<6; ++i){
    for (int j=0; j<8; ++j){
     a(i,j)= (double) ++c;
   }
  }

  auto b = pressio::ops::clone(a);
  ASSERT_EQ(b.rows(), 6);
  ASSERT_EQ(b.cols(), 8);
  for (int i=0; i<6; ++i){
    for (int j=0; j<8; ++j){
      ASSERT_DOUBLE_EQ(b(i,j),a(i,j));
   }
  }

  // check if b.data() == a.data()
  b(0, 0) = a(0, 0) + 1.;
  ASSERT_FALSE(b(0, 0) == a(0, 0));
}

TEST(ops_eigen_matrix, dense_matrix_extent)
{
  using T = Eigen::MatrixXd;
  T x(6,8);
  ASSERT_TRUE(pressio::ops::extent(x,0) == 6);
  ASSERT_TRUE(pressio::ops::extent(x,1) == 8);
  ASSERT_TRUE(pressio::ops::extent(x,2) == 1); // check extent over the rank
}

TEST(ops_eigen_matrix, dense_matrix_scale)
{
  using T = Eigen::MatrixXd;
  T a(6,8);
  a.setConstant(1.);

  pressio::ops::scale(a, 3.);
  for (int i=0; i<6; ++i){
   for (int j=0; j<8; ++j){
    ASSERT_DOUBLE_EQ(a(i,j),3.);
   }
  }
}

TEST(ops_eigen_matrix, dense_matrix_setzero)
{
  using T = Eigen::MatrixXd;
  T a(6,6);
  a.setConstant(11.);

  pressio::ops::set_zero(a);
  for (int i=0; i<6; ++i){
   for (int j=0; j<6; ++j){
    ASSERT_DOUBLE_EQ(a(i,j),0.);
   }
  }
}

TEST(ops_eigen_matrix, dense_matrix_fill)
{
  using T = Eigen::MatrixXd;
  T a(6,6);
  a.setConstant(11.);

  pressio::ops::fill(a, 55.);
  for (int i=0; i<6; ++i){
   for (int j=0; j<6; ++j){
    ASSERT_DOUBLE_EQ(a(i,j), 55.);
   }
  }
}

TEST(ops_eigen_matrix, dense_matrix_resize)
{
  using T = Eigen::MatrixXd;
  T a(6,6);
  pressio::ops::resize(a,3,4);
  ASSERT_EQ(a.rows(), 3);
  ASSERT_EQ(a.cols(), 4);
}

TEST(ops_eigen_matrix, dense_matrix_deep_copy)
{
  using T = Eigen::MatrixXd;
  T a(6,5);
  a.setConstant(44.);

  T b(6,5);
  pressio::ops::deep_copy(b,a);
  for (int i=0; i<6; ++i){
   for (int j=0; j<5; ++j){
    ASSERT_DOUBLE_EQ(b(i,j),44.);
   }
  }
}

TEST(ops_eigen_matrix, matrix_min_max)
{
  using T = Eigen::MatrixXd;
  T a(5, 5);
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      a(i, j) = 5 * i + j + 1;
    }
  }
  ASSERT_DOUBLE_EQ(pressio::ops::min(a), 1.);
  ASSERT_DOUBLE_EQ(pressio::ops::max(a), 25.);
}

TEST(ops_eigen_matrix, add_to_diagonal)
{
  using T = Eigen::MatrixXd;
  T A(6,6);
  A.setConstant(2.2);
  pressio::ops::add_to_diagonal(A, 1.1);
  for (int i=0; i<6; ++i){
    for (int j=0; j<6; ++j){
      if (i==j) {
        EXPECT_DOUBLE_EQ(A(i,j), 3.3);
      }
      else{
        EXPECT_DOUBLE_EQ(A(i,j), 2.2);
      }
    }
  }
}

TEST(ops_eigen_matrix, dense_matrix_update)
{
  Eigen::Matrix<double, 2, 2> M;
  Eigen::Matrix<double, 2, 2> A;
  M << 1., 2., 3., 4.;
  A << 5., 6., 7., 8.;

  pressio::ops::update(M, 2., A, 3.);
  EXPECT_DOUBLE_EQ(M(0, 0), 17.);
  EXPECT_DOUBLE_EQ(M(0, 1), 22.);
  EXPECT_DOUBLE_EQ(M(1, 0), 27.);
  EXPECT_DOUBLE_EQ(M(1, 1), 32.);

  // NaN injection through alpha=0
  const auto nan = std::nan("0");
  pressio::ops::fill(M, nan);
  pressio::ops::update(M, 0., A, 2.);
  EXPECT_DOUBLE_EQ(M(0, 0), 10.);
  EXPECT_DOUBLE_EQ(M(0, 1), 12.);
  EXPECT_DOUBLE_EQ(M(1, 0), 14.);
  EXPECT_DOUBLE_EQ(M(1, 1), 16.);

  // NaN injection through beta=0
  pressio::ops::fill(A, nan);
  pressio::ops::update(M, -1., A, 0.);
  EXPECT_DOUBLE_EQ(M(0, 0), -10.);
  EXPECT_DOUBLE_EQ(M(0, 1), -12.);
  EXPECT_DOUBLE_EQ(M(1, 0), -14.);
  EXPECT_DOUBLE_EQ(M(1, 1), -16.);

  // alpha=beta=0 corner case
  pressio::ops::fill(M, nan);
  pressio::ops::fill(A, nan);
  pressio::ops::update(M, 0., A, 0.);
  EXPECT_DOUBLE_EQ(M(0, 0), 0.);
  EXPECT_DOUBLE_EQ(M(0, 1), 0.);
  EXPECT_DOUBLE_EQ(M(1, 0), 0.);
  EXPECT_DOUBLE_EQ(M(1, 1), 0.);
}

TEST(ops_eigen_matrix, dense_matrix_update_epxr)
{
  Eigen::Matrix<double, 4, 4> M0;
  Eigen::Matrix<double, 4, 4> A0;
  pressio::ops::fill(M0, 1);
  pressio::ops::fill(A0, 2);
  auto M = pressio::subspan(M0, {1, 3}, {1, 3});
  auto A = pressio::subspan(A0, {1, 3}, {1, 3});

  pressio::ops::update(M, 2., A, 3.);
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      const bool sub = i > 0 && i < 3 && j > 0 && j < 3;
      EXPECT_DOUBLE_EQ(M0(i, j), sub ? 8.   // updated M
                                     : 1.); // unmodified part of M0
    }
  }
}
