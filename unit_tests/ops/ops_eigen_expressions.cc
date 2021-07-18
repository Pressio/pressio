
#include <gtest/gtest.h>
#include "pressio_ops.hpp"

TEST(ops_eigen, clone_vector)
{
  using T = Eigen::VectorXd;
  T a(6);
  for (int i=0; i<6; ++i){
   a(i)= (double) i;
  }

  auto b = pressio::ops::clone(a);
  ASSERT_EQ(b.size(), 6);
  ASSERT_FALSE( b.data()==a.data());
  for (int i=0; i<6; ++i){
    ASSERT_DOUBLE_EQ(b(i),a(i));
  }
}

TEST(ops_eigen, clone_matrix)
{
  using T = Eigen::MatrixXd;
  T a(6,8);

  int c=0;  
  for (int i=0; i<6; ++i){
    for (int j=0; j<8; ++j){
     a(i)= (double) ++c;
   }
  }

  auto b = pressio::ops::clone(a);
  ASSERT_EQ(b.rows(), 6);
  ASSERT_EQ(b.cols(), 8);
  ASSERT_FALSE( b.data()==a.data());

  for (int i=0; i<6; ++i){
    for (int j=0; j<8; ++j){
      ASSERT_DOUBLE_EQ(b(i,j),a(i,j));
     a(i)= (double) ++c;
   }
  }
}

TEST(ops_eigen, extent_vector)
{
  using vec_t = Eigen::VectorXd;
  vec_t x(6);
  ASSERT_TRUE(pressio::ops::extent(x,0)== 6);
}

TEST(ops_eigen, extent_matrix)
{
  using T = Eigen::MatrixXd;
  T x(6,8);
  ASSERT_TRUE(pressio::ops::extent(x,0) == 6);
  ASSERT_TRUE(pressio::ops::extent(x,1) == 8);
}

TEST(ops_eigen, abs)
{
  using vec_t = Eigen::VectorXd;
  vec_t x(6);
  for (int i=0; i<6; ++i){
    x(i) = -(double) i;
  }

  vec_t y(6);
  pressio::ops::abs(y,x);

  Eigen::VectorXd g(6);
  g(0) = 0.;
  g(1) = 1.;
  g(2) = 2.;
  g(3) = 3.;
  g(4) = 4.;
  g(5) = 5.;
  ASSERT_EQ( y.size(), 6 );
  for (int i=0; i<6; ++i){
    EXPECT_DOUBLE_EQ(y(i), g(i));
  }
}

TEST(ops_eigen, scale_vector)
{
  using T = Eigen::VectorXd;
  T a(6);
  a.setConstant(1.);

  pressio::ops::scale(a, 3.);
  for (int i=0; i<6; ++i){
    ASSERT_DOUBLE_EQ(a(i),3.);
  }
}

TEST(ops_eigen, scale_matrix)
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

TEST(ops_eigen, scale_span)
{
  using T = Eigen::VectorXd;
  T a(6);
  a.setConstant(1.);

  auto sp = pressio::span(a, 1,3);

  pressio::ops::scale(sp, 3.);
  ASSERT_DOUBLE_EQ(a(0),1.);
  ASSERT_DOUBLE_EQ(a(1),3.);
  ASSERT_DOUBLE_EQ(a(2),3.);
  ASSERT_DOUBLE_EQ(a(3),3.);
  ASSERT_DOUBLE_EQ(a(4),1.);
  ASSERT_DOUBLE_EQ(a(5),1.);
}


TEST(ops_eigen, scale_subspan)
{
  using T = Eigen::MatrixXd;
  T a(4,5);
  a.setConstant(1.);

  std::pair<int,int> r(1,3);
  std::pair<int,int> c(2,4);
  auto sp = pressio::subspan(a,r,c);

  pressio::ops::scale(sp, 3.);
  ASSERT_DOUBLE_EQ(a(0,0),1.);
  ASSERT_DOUBLE_EQ(a(0,1),1.);
  ASSERT_DOUBLE_EQ(a(0,2),1.);
  ASSERT_DOUBLE_EQ(a(0,3),1.);
  ASSERT_DOUBLE_EQ(a(0,4),1.);

  ASSERT_DOUBLE_EQ(a(1,0),1.);
  ASSERT_DOUBLE_EQ(a(1,1),1.);
  ASSERT_DOUBLE_EQ(a(1,2),3.);
  ASSERT_DOUBLE_EQ(a(1,3),3.);
  ASSERT_DOUBLE_EQ(a(1,4),1.);

  ASSERT_DOUBLE_EQ(a(2,0),1.);
  ASSERT_DOUBLE_EQ(a(2,1),1.);
  ASSERT_DOUBLE_EQ(a(2,2),3.);
  ASSERT_DOUBLE_EQ(a(2,3),3.);
  ASSERT_DOUBLE_EQ(a(2,4),1.);

  ASSERT_DOUBLE_EQ(a(3,0),1.);
  ASSERT_DOUBLE_EQ(a(3,1),1.);
  ASSERT_DOUBLE_EQ(a(3,2),1.);
  ASSERT_DOUBLE_EQ(a(3,3),1.);
  ASSERT_DOUBLE_EQ(a(3,4),1.);
}


TEST(ops_eigen, scale_diag)
{
  using T = Eigen::MatrixXd;
  T a(5,5);
  a.setConstant(1.);

  auto sp = pressio::diag(a);

  pressio::ops::scale(sp, 3.);
  ASSERT_DOUBLE_EQ(a(0,0),3.);
  ASSERT_DOUBLE_EQ(a(0,1),1.);
  ASSERT_DOUBLE_EQ(a(0,2),1.);
  ASSERT_DOUBLE_EQ(a(0,3),1.);
  ASSERT_DOUBLE_EQ(a(0,4),1.);

  ASSERT_DOUBLE_EQ(a(1,0),1.);
  ASSERT_DOUBLE_EQ(a(1,1),3.);
  ASSERT_DOUBLE_EQ(a(1,2),1.);
  ASSERT_DOUBLE_EQ(a(1,3),1.);
  ASSERT_DOUBLE_EQ(a(1,4),1.);

  ASSERT_DOUBLE_EQ(a(2,0),1.);
  ASSERT_DOUBLE_EQ(a(2,1),1.);
  ASSERT_DOUBLE_EQ(a(2,2),3.);
  ASSERT_DOUBLE_EQ(a(2,3),1.);
  ASSERT_DOUBLE_EQ(a(2,4),1.);

  ASSERT_DOUBLE_EQ(a(3,0),1.);
  ASSERT_DOUBLE_EQ(a(3,1),1.);
  ASSERT_DOUBLE_EQ(a(3,2),1.);
  ASSERT_DOUBLE_EQ(a(3,3),3.);
  ASSERT_DOUBLE_EQ(a(3,4),1.);

  ASSERT_DOUBLE_EQ(a(4,0),1.);
  ASSERT_DOUBLE_EQ(a(4,1),1.);
  ASSERT_DOUBLE_EQ(a(4,2),1.);
  ASSERT_DOUBLE_EQ(a(4,3),1.);
  ASSERT_DOUBLE_EQ(a(4,4),3.);
}

TEST(ops_eigen, setzero_vector)
{
  using T = Eigen::VectorXd;
  T a(6);
  a.setConstant(11.);

  pressio::ops::set_zero(a);
  for (int i=0; i<6; ++i){
    ASSERT_DOUBLE_EQ(a(i),0.);
  }
}


TEST(ops_eigen, setzero_matrix)
{
  using T = Eigen::VectorXd;
  T a(6);
  a.setConstant(11.);

  pressio::ops::set_zero(a);
  for (int i=0; i<6; ++i){
    ASSERT_DOUBLE_EQ(a(i),0.);
  }
}



