
#include <gtest/gtest.h>
#include "pressio/ops.hpp"

using vec_t = Eigen::VectorXd;
using mat_t = Eigen::MatrixXd;

TEST(ops_eigen, diag_extent)
{
  mat_t a(5,5);
  auto ex = pressio::diag(a);
  ASSERT_TRUE(pressio::ops::extent(ex,0)==5);
  ASSERT_TRUE(pressio::ops::extent(ex,1)==1); // check extent over the rank
}

TEST(ops_eigen, diag_abs)
{
  mat_t a(5,5);
  a.setConstant(-1);
  auto ex = pressio::diag(a);

  vec_t y(5);
  pressio::ops::abs(y,ex);
  ASSERT_DOUBLE_EQ(y(0),1.);
  ASSERT_DOUBLE_EQ(y(1),1.);
  ASSERT_DOUBLE_EQ(y(2),1.);
  ASSERT_DOUBLE_EQ(y(3),1.);
  ASSERT_DOUBLE_EQ(y(4),1.);
}

TEST(ops_eigen, diag_scale)
{
  using T = Eigen::MatrixXd;
  T a(5,5);
  a.setConstant(1.);

  auto exp = pressio::diag(a);

  pressio::ops::scale(exp, 3.);
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

TEST(ops_eigen, diag_set_zero)
{
  using T = Eigen::MatrixXd;
  T a(5,5);
  a.setConstant(1.);

  auto exp = pressio::diag(a);

  pressio::ops::set_zero(exp);
  ASSERT_DOUBLE_EQ(a(0,0),0.);
  ASSERT_DOUBLE_EQ(a(0,1),1.);
  ASSERT_DOUBLE_EQ(a(0,2),1.);
  ASSERT_DOUBLE_EQ(a(0,3),1.);
  ASSERT_DOUBLE_EQ(a(0,4),1.);

  ASSERT_DOUBLE_EQ(a(1,0),1.);
  ASSERT_DOUBLE_EQ(a(1,1),0.);
  ASSERT_DOUBLE_EQ(a(1,2),1.);
  ASSERT_DOUBLE_EQ(a(1,3),1.);
  ASSERT_DOUBLE_EQ(a(1,4),1.);

  ASSERT_DOUBLE_EQ(a(2,0),1.);
  ASSERT_DOUBLE_EQ(a(2,1),1.);
  ASSERT_DOUBLE_EQ(a(2,2),0.);
  ASSERT_DOUBLE_EQ(a(2,3),1.);
  ASSERT_DOUBLE_EQ(a(2,4),1.);

  ASSERT_DOUBLE_EQ(a(3,0),1.);
  ASSERT_DOUBLE_EQ(a(3,1),1.);
  ASSERT_DOUBLE_EQ(a(3,2),1.);
  ASSERT_DOUBLE_EQ(a(3,3),0.);
  ASSERT_DOUBLE_EQ(a(3,4),1.);

  ASSERT_DOUBLE_EQ(a(4,0),1.);
  ASSERT_DOUBLE_EQ(a(4,1),1.);
  ASSERT_DOUBLE_EQ(a(4,2),1.);
  ASSERT_DOUBLE_EQ(a(4,3),1.);
  ASSERT_DOUBLE_EQ(a(4,4),0.);
}

TEST(ops_eigen, diag_fill)
{
  using T = Eigen::MatrixXd;
  T a(5,5);
  a.setConstant(1.);

  auto exp = pressio::diag(a);

  pressio::ops::fill(exp, 44.);
  ASSERT_DOUBLE_EQ(a(0,0),44.);
  ASSERT_DOUBLE_EQ(a(0,1),1.);
  ASSERT_DOUBLE_EQ(a(0,2),1.);
  ASSERT_DOUBLE_EQ(a(0,3),1.);
  ASSERT_DOUBLE_EQ(a(0,4),1.);

  ASSERT_DOUBLE_EQ(a(1,0),1.);
  ASSERT_DOUBLE_EQ(a(1,1),44.);
  ASSERT_DOUBLE_EQ(a(1,2),1.);
  ASSERT_DOUBLE_EQ(a(1,3),1.);
  ASSERT_DOUBLE_EQ(a(1,4),1.);

  ASSERT_DOUBLE_EQ(a(2,0),1.);
  ASSERT_DOUBLE_EQ(a(2,1),1.);
  ASSERT_DOUBLE_EQ(a(2,2),44.);
  ASSERT_DOUBLE_EQ(a(2,3),1.);
  ASSERT_DOUBLE_EQ(a(2,4),1.);

  ASSERT_DOUBLE_EQ(a(3,0),1.);
  ASSERT_DOUBLE_EQ(a(3,1),1.);
  ASSERT_DOUBLE_EQ(a(3,2),1.);
  ASSERT_DOUBLE_EQ(a(3,3),44.);
  ASSERT_DOUBLE_EQ(a(3,4),1.);

  ASSERT_DOUBLE_EQ(a(4,0),1.);
  ASSERT_DOUBLE_EQ(a(4,1),1.);
  ASSERT_DOUBLE_EQ(a(4,2),1.);
  ASSERT_DOUBLE_EQ(a(4,3),1.);
  ASSERT_DOUBLE_EQ(a(4,4),44.);
}

TEST(ops_eigen, diag_min_max)
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
  auto exp = pressio::diag(a);
  ASSERT_DOUBLE_EQ(pressio::ops::min(exp), 0.);
  ASSERT_DOUBLE_EQ(pressio::ops::max(exp), 24.);
}

TEST(ops_eigen, diag_norms)
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
  gold << 0.,6.,12.,18.,24.;

  auto exp = pressio::diag(a);
  ASSERT_DOUBLE_EQ(pressio::ops::norm1(exp), gold.lpNorm<1>());
  ASSERT_DOUBLE_EQ(pressio::ops::norm2(exp), gold.lpNorm<2>());
}

TEST(ops_eigen, diag_dot_vector)
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

  auto exp = pressio::diag(a);
  ASSERT_DOUBLE_EQ(pressio::ops::dot(exp,b), 120.);
}

TEST(ops_eigen, diag_dot_diag)
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
  gold << 0.,6.,12.,18.,24.;

  auto exp = pressio::diag(a);
  ASSERT_DOUBLE_EQ(pressio::ops::dot(exp,exp),1080.);
}

TEST(ops_eigen, diag_pow)
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
  gold << 0.,6.,12.,18.,24.;

  auto exp = pressio::diag(a);
  pressio::ops::pow(exp, 2.);

  for (int i=0; i<a.rows(); ++i){
    EXPECT_DOUBLE_EQ(exp(i), gold(i)*gold(i));
  }
}

namespace {
Eigen::MatrixXd createMatrixForUpdate(){
  Eigen::MatrixXd M(3,3);
  M(0,0) = 1.;
  M(1,1) = 2.;
  M(2,2) = 3.;
  return M;
}
}

TEST(ops_eigen, diag_update1)
{
  auto M1 = createMatrixForUpdate();
  auto d1 = pressio::diag(M1);
  auto M2 = createMatrixForUpdate();
  auto d2 = pressio::diag(M2);

  pressio::ops::update(d1, 1., d2, 1.);
  EXPECT_DOUBLE_EQ( d1(0), 2.0);
  EXPECT_DOUBLE_EQ( d1(1), 4.0);
  EXPECT_DOUBLE_EQ( d1(2), 6.0);

  pressio::ops::update(d1, 0., d2, 1.);
  EXPECT_DOUBLE_EQ( d1(0), 1.0);
  EXPECT_DOUBLE_EQ( d1(1), 2.0);
  EXPECT_DOUBLE_EQ( d1(2), 3.0);
}

// TEST(ops_eigen, diag_update2)
// {
//   auto M1 = createMatrixForUpdate();
//   auto v = pressio::diag(M1);
//   auto M2 = createMatrixForUpdate();
//   auto a = pressio::diag(M2);
//   auto M3 = createMatrixForUpdate();
//   auto b = pressio::diag(M3);

//   pressio::ops::update(v, 1., a, 1., b, 1.);
//   EXPECT_DOUBLE_EQ( v(0), 3.0);
//   EXPECT_DOUBLE_EQ( v(1), 6.0);
//   EXPECT_DOUBLE_EQ( v(2), 9.0);

//   pressio::ops::update(v, 0., a, 1., b, 1.);
//   EXPECT_DOUBLE_EQ( v(0), 2.0);
//   EXPECT_DOUBLE_EQ( v(1), 4.0);
//   EXPECT_DOUBLE_EQ( v(2), 6.0);
// }

// TEST(ops_eigen, diag_update3)
// {
//   auto M1 = createMatrixForUpdate();
//   auto v = pressio::diag(M1);
//   auto M2 = createMatrixForUpdate();
//   auto a = pressio::diag(M2);
//   auto M3 = createMatrixForUpdate();
//   auto b = pressio::diag(M3);
//   auto M4 = createMatrixForUpdate();
//   auto c = pressio::diag(M4);

//   pressio::ops::update(v, 1., a, 1., b, 1., c, 1.);
//   EXPECT_DOUBLE_EQ( v(0), 4.0);
//   EXPECT_DOUBLE_EQ( v(1), 8.0);
//   EXPECT_DOUBLE_EQ( v(2), 12.0);

//   pressio::ops::update(v, 0., a, 1., b, 1., c, 1.);
//   EXPECT_DOUBLE_EQ( v(0), 3.0);
//   EXPECT_DOUBLE_EQ( v(1), 6.0);
//   EXPECT_DOUBLE_EQ( v(2), 9.0);
// }

// TEST(ops_eigen, diag_update4)
// {
//   auto M1 = createMatrixForUpdate();
//   auto v = pressio::diag(M1);
//   auto M2 = createMatrixForUpdate();
//   auto a = pressio::diag(M2);
//   auto M3 = createMatrixForUpdate();
//   auto b = pressio::diag(M3);
//   auto M4 = createMatrixForUpdate();
//   auto c = pressio::diag(M4);
//   auto M5 = createMatrixForUpdate();
//   auto d = pressio::diag(M5);

//   pressio::ops::update(v, 1., a, 1., b, 1., c, 1., d, 1.);
//   EXPECT_DOUBLE_EQ( v(0), 5.0);
//   EXPECT_DOUBLE_EQ( v(1), 10.0);
//   EXPECT_DOUBLE_EQ( v(2), 15.0);

//   pressio::ops::update(v, 0., a, 1., b, 1., c, 1., d, 1.);
//   EXPECT_DOUBLE_EQ( v(0), 4.0);
//   EXPECT_DOUBLE_EQ( v(1), 8.0);
//   EXPECT_DOUBLE_EQ( v(2), 12.0);
// }

TEST(ops_eigen, diag_elementwiseMultiply)
{
  Eigen::MatrixXd M1(3,3);
  M1(0,0)=1.; M1(1,1)=2.; M1(2,2)=3.;
  auto y = pressio::diag(M1);

  Eigen::MatrixXd M2(3,3);
  M2(0,0)=2.; M2(1,1)=3.; M2(2,2)=4.;
  const auto x = pressio::diag(M2);

  Eigen::MatrixXd M3(3,3);
  M3(0,0)=3.; M3(1,1)=4.; M3(2,2)=5.;
  const auto z = pressio::diag(M3);

  pressio::ops::elementwise_multiply(1., x, z, 1., y);
  EXPECT_DOUBLE_EQ( y(0), 7.0);
  EXPECT_DOUBLE_EQ( y(1), 14.0);
  EXPECT_DOUBLE_EQ( y(2), 23.0);
}
