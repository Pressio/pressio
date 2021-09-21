
#include <gtest/gtest.h>
#include "pressio/ops.hpp"


TEST(ops_eigen, subspan_extent)
{
  using T = Eigen::MatrixXd;
  T A(5,5);
  std::pair<int,int> r(1,3);
  std::pair<int,int> c(2,5);

  auto ex = pressio::subspan(A,r,c);
  ASSERT_TRUE(pressio::ops::extent(ex,0)==2);
  ASSERT_TRUE(pressio::ops::extent(ex,1)==3);
}

TEST(ops_eigen, subspan_scale)
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

TEST(ops_eigen, subspan_set_zero)
{
  using T = Eigen::MatrixXd;
  T a(4,5);
  a.setConstant(1.);

  std::pair<int,int> r(1,3);
  std::pair<int,int> c(2,4);
  auto sp = pressio::subspan(a,r,c);

  pressio::ops::set_zero(sp);
  ASSERT_DOUBLE_EQ(a(0,0),1.);
  ASSERT_DOUBLE_EQ(a(0,1),1.);
  ASSERT_DOUBLE_EQ(a(0,2),1.);
  ASSERT_DOUBLE_EQ(a(0,3),1.);
  ASSERT_DOUBLE_EQ(a(0,4),1.);

  ASSERT_DOUBLE_EQ(a(1,0),1.);
  ASSERT_DOUBLE_EQ(a(1,1),1.);
  ASSERT_DOUBLE_EQ(a(1,2),0.);
  ASSERT_DOUBLE_EQ(a(1,3),0.);
  ASSERT_DOUBLE_EQ(a(1,4),1.);

  ASSERT_DOUBLE_EQ(a(2,0),1.);
  ASSERT_DOUBLE_EQ(a(2,1),1.);
  ASSERT_DOUBLE_EQ(a(2,2),0.);
  ASSERT_DOUBLE_EQ(a(2,3),0.);
  ASSERT_DOUBLE_EQ(a(2,4),1.);

  ASSERT_DOUBLE_EQ(a(3,0),1.);
  ASSERT_DOUBLE_EQ(a(3,1),1.);
  ASSERT_DOUBLE_EQ(a(3,2),1.);
  ASSERT_DOUBLE_EQ(a(3,3),1.);
  ASSERT_DOUBLE_EQ(a(3,4),1.);
}

TEST(ops_eigen, subspan_fill)
{
  using T = Eigen::MatrixXd;
  T a(4,5);
  a.setConstant(1.);

  std::pair<int,int> r(1,3);
  std::pair<int,int> c(2,4);
  auto sp = pressio::subspan(a,r,c);

  pressio::ops::fill(sp, 44.);
  ASSERT_DOUBLE_EQ(a(0,0),1.);
  ASSERT_DOUBLE_EQ(a(0,1),1.);
  ASSERT_DOUBLE_EQ(a(0,2),1.);
  ASSERT_DOUBLE_EQ(a(0,3),1.);
  ASSERT_DOUBLE_EQ(a(0,4),1.);

  ASSERT_DOUBLE_EQ(a(1,0),1.);
  ASSERT_DOUBLE_EQ(a(1,1),1.);
  ASSERT_DOUBLE_EQ(a(1,2),44.);
  ASSERT_DOUBLE_EQ(a(1,3),44.);
  ASSERT_DOUBLE_EQ(a(1,4),1.);

  ASSERT_DOUBLE_EQ(a(2,0),1.);
  ASSERT_DOUBLE_EQ(a(2,1),1.);
  ASSERT_DOUBLE_EQ(a(2,2),44.);
  ASSERT_DOUBLE_EQ(a(2,3),44.);
  ASSERT_DOUBLE_EQ(a(2,4),1.);

  ASSERT_DOUBLE_EQ(a(3,0),1.);
  ASSERT_DOUBLE_EQ(a(3,1),1.);
  ASSERT_DOUBLE_EQ(a(3,2),1.);
  ASSERT_DOUBLE_EQ(a(3,3),1.);
  ASSERT_DOUBLE_EQ(a(3,4),1.);
}
