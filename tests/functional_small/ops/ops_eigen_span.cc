
#include <gtest/gtest.h>
#include "pressio/ops.hpp"

using vec_t = Eigen::VectorXd;
using mat_t = Eigen::MatrixXd;

TEST(ops_eigen, span_extent)
{
  vec_t a(8);
  auto ex = pressio::span(a,5,2);
  ASSERT_TRUE(pressio::ops::extent(ex,0)==2);
}

TEST(ops_eigen, span_abs)
{
  vec_t a(8);
  a.setConstant(-1);
  auto ex = pressio::span(a,5,2);

  vec_t y(5);
  pressio::ops::abs(y,ex);
  ASSERT_DOUBLE_EQ(y(0),1.);
  ASSERT_DOUBLE_EQ(y(1),1.);
}

TEST(ops_eigen, span_scale)
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

TEST(ops_eigen, span_set_zero)
{
  using T = Eigen::VectorXd;
  T a(6);
  a.setConstant(1.2);
  auto sp = pressio::span(a, 1,3);
  pressio::ops::set_zero(sp);
  ASSERT_DOUBLE_EQ(a(0),1.2);
  ASSERT_DOUBLE_EQ(a(1),0.);
  ASSERT_DOUBLE_EQ(a(2),0.);
  ASSERT_DOUBLE_EQ(a(3),0.);
  ASSERT_DOUBLE_EQ(a(4),1.2);
  ASSERT_DOUBLE_EQ(a(5),1.2);
}

TEST(ops_eigen, span_fill)
{
  using T = Eigen::VectorXd;
  T a(6);
  a.setConstant(1.2);
  auto sp = pressio::span(a, 1,3);
  pressio::ops::fill(sp, 44.);
  ASSERT_DOUBLE_EQ(a(0),1.2);
  ASSERT_DOUBLE_EQ(a(1),44.);
  ASSERT_DOUBLE_EQ(a(2),44.);
  ASSERT_DOUBLE_EQ(a(3),44.);
  ASSERT_DOUBLE_EQ(a(4),1.2);
  ASSERT_DOUBLE_EQ(a(5),1.2);
}

TEST(ops_eigen, span_min_max)
{
  using T = Eigen::VectorXd;
  T a(6);
  for (int i=0; i<6; ++i){
    a(i) = (double) i;
  }
  auto sp = pressio::span(a, 1,3);
  ASSERT_DOUBLE_EQ(pressio::ops::min(sp), 1.);
  ASSERT_DOUBLE_EQ(pressio::ops::max(sp), 3.);
}

TEST(ops_eigen, span_norms)
{
  using T = Eigen::VectorXd;
  T a(6);
  for (int i=0; i<6; ++i){
    a(i) = (double) i;
  }
  auto sp = pressio::span(a, 1,3);
  ASSERT_DOUBLE_EQ(pressio::ops::norm1(sp), 6.);
  ASSERT_DOUBLE_EQ(pressio::ops::norm2(sp), std::sqrt(14.));
}

TEST(ops_eigen, span_dot)
{
  using T = Eigen::VectorXd;
  T a(6);
  for (int i=0; i<6; ++i){
    a(i) = (double) i;
  }
  auto sp = pressio::span(a, 1,3);

  Eigen::VectorXd b(3);
  b.setConstant(2.);

  ASSERT_DOUBLE_EQ(pressio::ops::dot(sp, b), 12.);
}

TEST(ops_eigen, span_pow)
{
  using vec_t = Eigen::VectorXd;
  vec_t x(8);
  for (int i=2; i<8; ++i) x(i) = (double) (i-2);

  auto sp = pressio::span(x, 2,6);
  pressio::ops::pow(sp, 2.);
  Eigen::VectorXd g(6);
  g(0) = 0.;
  g(1) = 1.;
  g(2) = 4.;
  g(3) = 9.;
  g(4) = 16.;
  g(5) = 25.;
  for (int i=0; i<6; ++i){
    EXPECT_DOUBLE_EQ(sp(i), g(i));
  }
}


namespace {
Eigen::VectorXd createVectorForUpdate(){
  Eigen::VectorXd M(6);
  M(2) = 1.;
  M(3) = 2.;
  M(4) = 3.;
  return M;
}
}

TEST(ops_eigen, span_update1)
{
  auto M1 = createVectorForUpdate();
  auto d1 = pressio::span(M1,2,3);
  auto M2 = createVectorForUpdate();
  auto d2 = pressio::span(M2,2,3);

  pressio::ops::update(d1, 1., d2, 1.);
  EXPECT_DOUBLE_EQ( d1(0), 2.0);
  EXPECT_DOUBLE_EQ( d1(1), 4.0);
  EXPECT_DOUBLE_EQ( d1(2), 6.0);

  pressio::ops::update(d1, 0., d2, 1.);
  EXPECT_DOUBLE_EQ( d1(0), 1.0);
  EXPECT_DOUBLE_EQ( d1(1), 2.0);
  EXPECT_DOUBLE_EQ( d1(2), 3.0);
}

TEST(ops_eigen, span_update2)
{
  auto M1 = createVectorForUpdate();
  auto v = pressio::span(M1,2,3);
  auto M2 = createVectorForUpdate();
  auto a = pressio::span(M2,2,3);
  auto M3 = createVectorForUpdate();
  auto b = pressio::span(M3,2,3);

  pressio::ops::update(v, 1., a, 1., b, 1.);
  EXPECT_DOUBLE_EQ( v(0), 3.0);
  EXPECT_DOUBLE_EQ( v(1), 6.0);
  EXPECT_DOUBLE_EQ( v(2), 9.0);

  pressio::ops::update(v, 0., a, 1., b, 1.);
  EXPECT_DOUBLE_EQ( v(0), 2.0);
  EXPECT_DOUBLE_EQ( v(1), 4.0);
  EXPECT_DOUBLE_EQ( v(2), 6.0);
}

TEST(ops_eigen, span_update3)
{
  auto M1 = createVectorForUpdate();
  auto v = pressio::span(M1,2,3);
  auto M2 = createVectorForUpdate();
  auto a = pressio::span(M2,2,3);
  auto M3 = createVectorForUpdate();
  auto b = pressio::span(M3,2,3);
  auto M4 = createVectorForUpdate();
  auto c = pressio::span(M4,2,3);

  pressio::ops::update(v, 1., a, 1., b, 1., c, 1.);
  EXPECT_DOUBLE_EQ( v(0), 4.0);
  EXPECT_DOUBLE_EQ( v(1), 8.0);
  EXPECT_DOUBLE_EQ( v(2), 12.0);

  pressio::ops::update(v, 0., a, 1., b, 1., c, 1.);
  EXPECT_DOUBLE_EQ( v(0), 3.0);
  EXPECT_DOUBLE_EQ( v(1), 6.0);
  EXPECT_DOUBLE_EQ( v(2), 9.0);
}

TEST(ops_eigen, span_update4)
{
  auto M1 = createVectorForUpdate();
  auto v = pressio::span(M1,2,3);
  auto M2 = createVectorForUpdate();
  auto a = pressio::span(M2,2,3);
  auto M3 = createVectorForUpdate();
  auto b = pressio::span(M3,2,3);
  auto M4 = createVectorForUpdate();
  auto c = pressio::span(M4,2,3);
  auto M5 = createVectorForUpdate();
  auto d = pressio::span(M5,2,3);

  pressio::ops::update(v, 1., a, 1., b, 1., c, 1., d, 1.);
  EXPECT_DOUBLE_EQ( v(0), 5.0);
  EXPECT_DOUBLE_EQ( v(1), 10.0);
  EXPECT_DOUBLE_EQ( v(2), 15.0);

  pressio::ops::update(v, 0., a, 1., b, 1., c, 1., d, 1.);
  EXPECT_DOUBLE_EQ( v(0), 4.0);
  EXPECT_DOUBLE_EQ( v(1), 8.0);
  EXPECT_DOUBLE_EQ( v(2), 12.0);
}

TEST(ops_eigen, span_elementwiseMultiply)
{
  Eigen::VectorXd M1(6);
  M1(2)=1.; M1(3)=2.; M1(4)=3.;
  auto y = pressio::span(M1,2,3);

  Eigen::VectorXd M2(6);
  M2(2)=2.; M2(3)=3.; M2(4)=4.;
  const auto x = pressio::span(M2,2,3);

  Eigen::VectorXd M3(6);
  M3(2)=3.; M3(3)=4.; M3(4)=5.;
  const auto z = pressio::span(M3,2,3);

  pressio::ops::elementwise_multiply(1., x, z, 1., y);
  EXPECT_DOUBLE_EQ( y(0), 7.0);
  EXPECT_DOUBLE_EQ( y(1), 14.0);
  EXPECT_DOUBLE_EQ( y(2), 23.0);
}
