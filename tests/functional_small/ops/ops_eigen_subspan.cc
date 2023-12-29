
#include <gtest/gtest.h>
#include "pressio/ops.hpp"


TEST(ops_eigen_subspan, extent)
{
  using T = Eigen::MatrixXd;
  T A(5,5);
  std::pair<int,int> r(1,3);
  std::pair<int,int> c(2,5);

  auto ex = pressio::subspan(A,r,c);
  ASSERT_EQ(2, pressio::ops::extent(ex, 0));
  ASSERT_EQ(3, pressio::ops::extent(ex, 1));
  ASSERT_EQ(1, pressio::ops::extent(ex, 2)); // check extent over the rank
}

TEST(ops_eigen_subspan, scale)
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

TEST(ops_eigen_subspan, set_zero)
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

TEST(ops_eigen_subspan, fill)
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

TEST(ops_eigen_subspan, min_max)
{
  using T = Eigen::MatrixXd;
  T a(5, 5);
  for (int i=0; i<5; ++i){
    for (int j=0; j<5; ++j){
      a(i, j)= (double)(i * 5 + j);
    }
  }

  std::pair<int,int> r(1,3);
  std::pair<int,int> c(2,4);
  auto sp = pressio::subspan(a,r,c);

  ASSERT_DOUBLE_EQ(pressio::ops::min(sp), 7.);
  ASSERT_DOUBLE_EQ(pressio::ops::max(sp), 13.);
}
