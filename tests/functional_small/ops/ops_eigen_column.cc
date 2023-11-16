
#include <gtest/gtest.h>
#include "pressio/ops.hpp"

using vec_t = Eigen::VectorXd;
using mat_t = Eigen::MatrixXd;

TEST(ops_eigen, column_extent)
{
  mat_t a(5,5);
  auto ex = pressio::column(a, 0);
  ASSERT_TRUE(pressio::ops::extent(ex,0)==5);
  ASSERT_TRUE(pressio::ops::extent(ex,1)==1);
}

TEST(ops_eigen, column_abs)
{
  mat_t a(5,5);
  a.setConstant(-1);
  auto ex = pressio::column(a, 0);

  vec_t y(5);
  pressio::ops::abs(y,ex);
  ASSERT_DOUBLE_EQ(y(0),1.);
  ASSERT_DOUBLE_EQ(y(1),1.);
  ASSERT_DOUBLE_EQ(y(2),1.);
  ASSERT_DOUBLE_EQ(y(3),1.);
  ASSERT_DOUBLE_EQ(y(4),1.);
}

TEST(ops_eigen, column_scale)
{
  using T = Eigen::MatrixXd;
  T a(5,5);
  a.setConstant(1.);

  auto exp = pressio::column(a, 1);
  pressio::ops::scale(exp, 3.);

  for (auto j : {0,2,3,4}){
    ASSERT_DOUBLE_EQ(a(0,j),1.);
    ASSERT_DOUBLE_EQ(a(1,j),1.);
    ASSERT_DOUBLE_EQ(a(2,j),1.);
    ASSERT_DOUBLE_EQ(a(3,j),1.);
    ASSERT_DOUBLE_EQ(a(4,j),1.);
  }

  ASSERT_DOUBLE_EQ(a(0,1),3.);
  ASSERT_DOUBLE_EQ(a(1,1),3.);
  ASSERT_DOUBLE_EQ(a(2,1),3.);
  ASSERT_DOUBLE_EQ(a(3,1),3.);
  ASSERT_DOUBLE_EQ(a(4,1),3.);
}

TEST(ops_eigen, column_set_zero)
{
  using T = Eigen::MatrixXd;
  T a(5,5);
  a.setConstant(1.);

  auto exp = pressio::column(a, 1);
  pressio::ops::set_zero(exp);

  for (auto j : {0,2,3,4}){
    ASSERT_DOUBLE_EQ(a(0,j),1.);
    ASSERT_DOUBLE_EQ(a(1,j),1.);
    ASSERT_DOUBLE_EQ(a(2,j),1.);
    ASSERT_DOUBLE_EQ(a(3,j),1.);
    ASSERT_DOUBLE_EQ(a(4,j),1.);
  }

  ASSERT_DOUBLE_EQ(a(0,1),0.);
  ASSERT_DOUBLE_EQ(a(1,1),0.);
  ASSERT_DOUBLE_EQ(a(2,1),0.);
  ASSERT_DOUBLE_EQ(a(3,1),0.);
  ASSERT_DOUBLE_EQ(a(4,1),0.);
}

TEST(ops_eigen, column_fill)
{
  using T = Eigen::MatrixXd;
  T a(5,5);
  a.setConstant(1.);

  auto exp = pressio::column(a, 1);
  pressio::ops::fill(exp, 44.);

  for (auto j : {0,2,3,4}){
    ASSERT_DOUBLE_EQ(a(0,j),1.);
    ASSERT_DOUBLE_EQ(a(1,j),1.);
    ASSERT_DOUBLE_EQ(a(2,j),1.);
    ASSERT_DOUBLE_EQ(a(3,j),1.);
    ASSERT_DOUBLE_EQ(a(4,j),1.);
  }

  ASSERT_DOUBLE_EQ(a(0,1),44.);
  ASSERT_DOUBLE_EQ(a(1,1),44.);
  ASSERT_DOUBLE_EQ(a(2,1),44.);
  ASSERT_DOUBLE_EQ(a(3,1),44.);
  ASSERT_DOUBLE_EQ(a(4,1),44.);
}

TEST(ops_eigen, column_min_max)
{
  using T = Eigen::MatrixXd;
  T a(5,5);

  double v = 0;
  for (int i=0; i<5; ++i){
   for (int j=0; j<5; ++j){
    a(i,j) = v;
    v+=1.;
   }
  }

  auto exp = pressio::column(a, 1);
  ASSERT_DOUBLE_EQ(pressio::ops::min(exp), 1.);
  ASSERT_DOUBLE_EQ(pressio::ops::max(exp), 21.);
}

TEST(ops_eigen, column_norms)
{
  using T = Eigen::MatrixXd;
  T a(5,5);

  double v = 0;
  for (int i=0; i<5; ++i){
   for (int j=0; j<5; ++j){
    a(i,j) = v;
    v+=1.;
   }
  }

  Eigen::VectorXd gold(5);
  gold << 1.,6.,11.,16.,21.;

  auto exp = pressio::column(a, 1);
  ASSERT_DOUBLE_EQ(pressio::ops::norm1(exp), gold.lpNorm<1>());
  ASSERT_DOUBLE_EQ(pressio::ops::norm2(exp), gold.lpNorm<2>());
}

TEST(ops_eigen, column_dot_vector)
{
  using T = Eigen::MatrixXd;
  T a(5,5);

  double v = 0;
  for (int i=0; i<5; ++i){
   for (int j=0; j<5; ++j){
    a(i,j) = v;
    v+=1.;
   }
  }

  Eigen::VectorXd b(5);
  b.setConstant(2.);

  auto exp = pressio::column(a, 1);
  ASSERT_DOUBLE_EQ(pressio::ops::dot(exp,b), 110.);
}

TEST(ops_eigen, column_dot_diag)
{
  using T = Eigen::MatrixXd;
  T a(5,5);
  // 0,1,2,3,4
  // 5,6,7,8,9
  // 10,11,12,13,14
  // 15,16,17,18,19
  // 20,21,22,23,24
  double v = 0;
  for (int i=0; i<5; ++i){
   for (int j=0; j<5; ++j){
    a(i,j) = v;
    v+=1.;
   }
  }

  Eigen::VectorXd gold(5);
  gold << 1.,6.,11.,16.,21.;

  auto exp = pressio::column(a,1);
  ASSERT_DOUBLE_EQ(pressio::ops::dot(exp,exp),855.);
}

TEST(ops_eigen, column_pow)
{
  using T = Eigen::MatrixXd;
  T a(5,5);
  double v = 0;
  for (int i=0; i<5; ++i){
   for (int j=0; j<5; ++j){
    a(i,j) = v;
    v+=1.;
   }
  }

  Eigen::VectorXd gold(5);
  gold << 1.,6.,11.,16.,21.;

  auto exp = pressio::column(a,1);
  pressio::ops::pow(exp, 2.);

  for (int i=0; i<a.rows(); ++i){
    EXPECT_DOUBLE_EQ(exp(i), gold(i)*gold(i));
  }
}

namespace {
Eigen::MatrixXd createMatrixForUpdate(){
  Eigen::MatrixXd M(3,3);
  M.col(0).setConstant(1.);
  M.col(1).setConstant(2.);
  M.col(2).setConstant(3.);
  return M;
}
}

TEST(ops_eigen, column_update1)
{
  auto M1 = createMatrixForUpdate();
  auto d1 = pressio::column(M1,0);
  auto M2 = createMatrixForUpdate();
  auto d2 = pressio::column(M2,1);

  pressio::ops::update(d1, 1., d2, 1.);
  EXPECT_DOUBLE_EQ( d1(0), 3.0);
  EXPECT_DOUBLE_EQ( d1(1), 3.0);
  EXPECT_DOUBLE_EQ( d1(2), 3.0);

  pressio::ops::update(d1, 0., d2, 1.);
  EXPECT_DOUBLE_EQ( d1(0), 2.0);
  EXPECT_DOUBLE_EQ( d1(1), 2.0);
  EXPECT_DOUBLE_EQ( d1(2), 2.0);
}
