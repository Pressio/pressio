
#include <gtest/gtest.h>
#include "pressio/ops.hpp"

TEST(ops_eigen, dense_matrix_clone)
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
  ASSERT_FALSE( b.data()==a.data());

  for (int i=0; i<6; ++i){
    for (int j=0; j<8; ++j){
      ASSERT_DOUBLE_EQ(b(i,j),a(i,j));
   }
  }
}

TEST(ops_eigen, dense_matrix_extent)
{
  using T = Eigen::MatrixXd;
  T x(6,8);
  ASSERT_TRUE(pressio::ops::extent(x,0) == 6);
  ASSERT_TRUE(pressio::ops::extent(x,1) == 8);
}

TEST(ops_eigen, dense_matrix_scale)
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

TEST(ops_eigen, dense_matrix_setzero)
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

TEST(ops_eigen, dense_matrix_fill)
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

TEST(ops_eigen, dense_matrix_resize)
{
  using T = Eigen::MatrixXd;
  T a(6,6);
  pressio::ops::resize(a,3,4);
  ASSERT_EQ(a.rows(), 3);
  ASSERT_EQ(a.cols(), 4);
}

TEST(ops_eigen, dense_matrix_deep_copy)
{
  using T = Eigen::MatrixXd;
  T a(6,5);
  a.setConstant(44.);

  T b(6,5);
  pressio::ops::deep_copy(b,a);
  for (int i=0; i<6; ++i){
   for (int j=0; j<5; ++j){
    ASSERT_DOUBLE_EQ(a(i,j),44.);
   }
  }
}
