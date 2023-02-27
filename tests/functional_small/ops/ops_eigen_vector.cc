
#include <gtest/gtest.h>
#include "pressio/ops.hpp"

using T = Eigen::VectorXd;
using V_t = Eigen::Matrix<double,3,1>;

TEST(ops_eigen, vector_clone)
{
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

TEST(ops_eigen, vector_extent)
{
  T x(6);
  ASSERT_TRUE(pressio::ops::extent(x,0)== 6);
}

TEST(ops_eigen, vector_abs)
{
  T x(6);
  for (int i=0; i<6; ++i){
    x(i) = -(double) i;
  }

  T y(6);
  pressio::ops::abs(y,x);

  T g(6);
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

TEST(ops_eigen, vector_scale)
{
  const size_t n = 6;
  T a(n);
  a.setConstant(2.);

  pressio::ops::scale(a, 3.);
  for (int i = 0; i < n; ++i){
    ASSERT_DOUBLE_EQ(a(i), 6.);
  }

  a.setConstant(std::nan("0"));
  pressio::ops::scale(a, 0.);
  for (int i = 0; i < n; ++i){
    ASSERT_DOUBLE_EQ(a(i), 0.);
  }
}

TEST(ops_eigen, vector_setzero)
{
  T a(6);
  a.setConstant(11.);

  pressio::ops::set_zero(a);
  for (int i=0; i<6; ++i){
    ASSERT_DOUBLE_EQ(a(i),0.);
  }
}

TEST(ops_eigen, vector_fill)
{
  T a(6);
  a.setConstant(11.);

  pressio::ops::fill(a, 44.88);
  for (int i=0; i<6; ++i){
    ASSERT_DOUBLE_EQ(a(i),44.88);
  }
}

TEST(ops_eigen, vector_resize)
{
  T a(6);
  pressio::ops::resize(a,3);
  ASSERT_EQ(a.size(), 3);
}

TEST(ops_eigen, vector_deep_copy)
{
  T a(6);
  a.setConstant(44.);

  T b(6);
  pressio::ops::deep_copy(b,a);
  for (int i=0; i<6; ++i){
    ASSERT_DOUBLE_EQ(a(i),44.);
  }
}

TEST(ops_eigen, vector_min_max)
{
  T a(5);
  for (int i=0; i<5; ++i){
   a(i)= (double) i;
  }

  ASSERT_DOUBLE_EQ(pressio::ops::min(a), 0.);
  ASSERT_DOUBLE_EQ(pressio::ops::max(a), 4.);
}

TEST(ops_eigen, vector_norm1)
{
  T a(5);
  for (int i=0; i<5; ++i){
   a(i)= (double) i;
  }
  std::cout << pressio::ops::norm1(a) << std::endl;
  ASSERT_DOUBLE_EQ(pressio::ops::norm1(a), a.lpNorm<1>());
}

TEST(ops_eigen, vector_norm2)
{
  T a(5);
  for (int i=0; i<5; ++i){
   a(i)= (double) i;
  }

  std::cout << pressio::ops::norm2(a) << std::endl;

 ASSERT_DOUBLE_EQ(pressio::ops::norm2(a), a.norm());
 ASSERT_DOUBLE_EQ(pressio::ops::norm2(a), a.lpNorm<2>());
}

TEST(ops_eigen, vector_dot)
{
  T a(5);
  a.setConstant(1.);
  T b(5);
  b.setConstant(2.);

  ASSERT_DOUBLE_EQ(pressio::ops::dot(a,b), 10.);
}

TEST(ops_eigen, vector_pow)
{
  T x(6);
  for (int i=0; i<6; ++i) x(i) = (double) i;

  ::pressio::ops::pow(x, 2.);
  Eigen::VectorXd g(6);
  g(0) = 0.;
  g(1) = 1.;
  g(2) = 4.;
  g(3) = 9.;
  g(4) = 16.;
  g(5) = 25.;
  ASSERT_EQ( x.size(), 6 );
  for (int i=0; i<6; ++i)
    EXPECT_DOUBLE_EQ(x(i), g(i));
}

TEST(ops_eigen, vector_absPowPos)
{
  T y(6);
  T x(6);
  for (int i=0; i<6; ++i){
    x(i) = (double) i; x(i)*=-1.;
  }

  ::pressio::ops::abs_pow(y, x, 3.);
  T g(6);
  g(0) = 0.;
  g(1) = 1.;
  g(2) = 8.;
  g(3) = 27.;
  g(4) = 64.;
  g(5) = 125.;
  ASSERT_EQ( y.size(), 6 );
  for (int i=0; i<6; ++i){
    EXPECT_DOUBLE_EQ(y(i), g(i));
  }
}

TEST(ops_eigen, vector_absPowNeg)
{
  T y(6);
  T x(6);
  for (int i=0; i<6; ++i){
    x(i) = (double) i;
    x(i)*=-1.;
  }

  ::pressio::ops::abs_pow(y, x, -3., 0.00001);
  // std::cout << y << std::endl;

  T g(6);
  g(0) = 1./0.00001; // because we guard against diving by zero
  g(1) = 1.;
  g(2) = 1./8.;
  g(3) = 1./27.;
  g(4) = 1./64.;
  g(5) = 1./125.;
  ASSERT_EQ( y.size(), 6 );
  for (int i=0; i<6; ++i){
    EXPECT_DOUBLE_EQ(y(i), g(i));
  }
}

TEST(ops_eigen, vector_update1)
{
  V_t v; v << 1.,2.,3.;
  V_t a; a << 1.,2.,3.;

  pressio::ops::update(v, 1., a, 1.);
  EXPECT_DOUBLE_EQ( v(0), 2.0);
  EXPECT_DOUBLE_EQ( v(1), 4.0);
  EXPECT_DOUBLE_EQ( v(2), 6.0);

  pressio::ops::update(v, 0., a, 1.);
  EXPECT_DOUBLE_EQ( v(0), 1.0);
  EXPECT_DOUBLE_EQ( v(1), 2.0);
  EXPECT_DOUBLE_EQ( v(2), 3.0);
}

TEST(ops_eigen, vector_update2)
{
  V_t v; v << 1.,2.,3.;
  V_t a; a << 1.,2.,3.;
  V_t b; b << 1.,2.,3.;

  pressio::ops::update(v, 1., a, 1., b, 1.);
  EXPECT_DOUBLE_EQ( v(0), 3.0);
  EXPECT_DOUBLE_EQ( v(1), 6.0);
  EXPECT_DOUBLE_EQ( v(2), 9.0);

  pressio::ops::update(v, 0., a, 1., b, 1.);
  EXPECT_DOUBLE_EQ( v(0), 2.0);
  EXPECT_DOUBLE_EQ( v(1), 4.0);
  EXPECT_DOUBLE_EQ( v(2), 6.0);
}

TEST(ops_eigen, vector_update3)
{
  V_t v; v << 1.,2.,3.;
  V_t a; a << 1.,2.,3.;
  V_t b; b << 1.,2.,3.;
  V_t c; c << 1.,2.,3.;

  pressio::ops::update(v, 1., a, 1., b, 1., c, 1.);
  EXPECT_DOUBLE_EQ( v(0), 4.0);
  EXPECT_DOUBLE_EQ( v(1), 8.0);
  EXPECT_DOUBLE_EQ( v(2), 12.0);

  pressio::ops::update(v, 0., a, 1., b, 1., c, 1.);
  EXPECT_DOUBLE_EQ( v(0), 3.0);
  EXPECT_DOUBLE_EQ( v(1), 6.0);
  EXPECT_DOUBLE_EQ( v(2), 9.0);
}

TEST(ops_eigen, vector_update4)
{
  V_t v; v << 1.,2.,3.;
  V_t a; a << 1.,2.,3.;
  V_t b; b << 1.,2.,3.;
  V_t c; c << 1.,2.,3.;
  V_t d; d << 1.,2.,3.;

  pressio::ops::update(v, 1., a, 1., b, 1., c, 1., d, 1.);
  EXPECT_DOUBLE_EQ( v(0), 5.0);
  EXPECT_DOUBLE_EQ( v(1), 10.0);
  EXPECT_DOUBLE_EQ( v(2), 15.0);

  pressio::ops::update(v, 0., a, 1., b, 1., c, 1., d, 1.);
  EXPECT_DOUBLE_EQ( v(0), 4.0);
  EXPECT_DOUBLE_EQ( v(1), 8.0);
  EXPECT_DOUBLE_EQ( v(2), 12.0);
}

TEST(ops_eigen, vector_update_nan1)
{
  V_t v, a, vecOfNans;
  pressio::ops::fill(v, 1.);
  pressio::ops::fill(a, 1.);
  pressio::ops::fill(vecOfNans, std::nan("0"));

  // Note: this test covers just enough nan/non-nan combinations
  // to trigger and verify all execution paths in our update()
  // implementations, which include anti-NaN-injection variants
  pressio::ops::update(v, 1., vecOfNans, 0.);
  EXPECT_DOUBLE_EQ(v(0), 1.0);

  pressio::ops::update(v, 1., vecOfNans, 0., vecOfNans, 0.);
  EXPECT_DOUBLE_EQ(v(0), 1.0);
  pressio::ops::update(v, 1., a, 1., vecOfNans, 0.);
  EXPECT_DOUBLE_EQ(v(0), 2.);

  pressio::ops::update(v, 1., vecOfNans, 0., vecOfNans, 0., vecOfNans, 0.);
  EXPECT_DOUBLE_EQ(v(0), 2.0);
  pressio::ops::update(v, 1., a, 1., vecOfNans, 0., a, 1.);
  EXPECT_DOUBLE_EQ(v(0), 4.);
  pressio::ops::update(v, 1., a, 1., a, 1., vecOfNans, 0.);
  EXPECT_DOUBLE_EQ(v(0), 6.);

  pressio::ops::update(v, 1., vecOfNans, 0., vecOfNans, 0., vecOfNans, 0., vecOfNans, 0.);
  EXPECT_DOUBLE_EQ(v(0), 6.0);
  pressio::ops::update(v, 1., a, 1., vecOfNans, 0., a, 1., a, 1.);
  EXPECT_DOUBLE_EQ(v(0), 9.);
  pressio::ops::update(v, 1., a, 1., a, 1., vecOfNans, 0., a, 1.);
  EXPECT_DOUBLE_EQ(v(0), 12.);
  pressio::ops::update(v, 1., a, 1., a, 1., a, 1., vecOfNans, 0.);
  EXPECT_DOUBLE_EQ(v(0), 15.);
}

// injects NaN through the updated vector
TEST(ops_eigen, vector_update_nan2)
{
  const auto nan = std::nan("0");
  V_t v, a;
  pressio::ops::fill(a, 1.);

  pressio::ops::fill(v, nan);
  pressio::ops::update(v, 0., a, 1.);
  EXPECT_DOUBLE_EQ(v(0), 1.0);

  pressio::ops::fill(v, nan);
  pressio::ops::update(v, 0., a, 0.);
  EXPECT_DOUBLE_EQ(v(0), 0.0);

  pressio::ops::fill(v, nan);
  pressio::ops::update(v, 0., a, 1., a, 1.);
  EXPECT_DOUBLE_EQ(v(0), 2.0);

  pressio::ops::fill(v, nan);
  pressio::ops::update(v, 0., a, 1., a, 1., a, 1.);
  EXPECT_DOUBLE_EQ(v(0), 3.0);

  pressio::ops::fill(v, nan);
  pressio::ops::update(v, 0., a, 1., a, 1., a, 1., a, 1.);
  EXPECT_DOUBLE_EQ(v(0), 4.0);
}

TEST(ops_eigen, vector_update_expr_span)
{
  T v0(5);
  T a0(5);
  auto v = pressio::span(v0, 1, 3);
  auto a = pressio::span(a0, 1, 3);
  pressio::ops::fill(v, 10.);
  pressio::ops::fill(a, 1.);

  pressio::ops::update(v, 1., a, 1.);
  pressio::ops::update(v, 1., a, 1., a, 2.);
  pressio::ops::update(v, 1., a, 1., a, 2., a, 3.);
  pressio::ops::update(v, 1., a, 1., a, 2., a, 3., a, 4.);

  // Note: just check the final result as this test is more about
  //       whether expressions compile and work than computation itself
  EXPECT_DOUBLE_EQ(v(0), 30.0);
  EXPECT_DOUBLE_EQ(v(1), 30.0);
  EXPECT_DOUBLE_EQ(v(2), 30.0);
}

TEST(ops_eigen, vector_update_expr_diag)
{
  Eigen::Matrix<double, 3, 3> v0;
  Eigen::Matrix<double, 3, 3> a0;
  auto v = pressio::diag(v0);
  auto a = pressio::diag(a0);
  pressio::ops::fill(v, 10.);
  pressio::ops::fill(a, 1.);

  pressio::ops::update(v, 1., a, 1.);
  pressio::ops::update(v, 1., a, 1., a, 2.);
  pressio::ops::update(v, 1., a, 1., a, 2., a, 3.);
  pressio::ops::update(v, 1., a, 1., a, 2., a, 3., a, 4.);

  // Note: just check the final result as this test is more about
  //       whether expressions compile and work than computation itself
  EXPECT_DOUBLE_EQ(v(0), 30.0);
  EXPECT_DOUBLE_EQ(v(1), 30.0);
  EXPECT_DOUBLE_EQ(v(2), 30.0);
}

TEST(ops_eigen, vector_elementwiseMultiply)
{
  V_t y; y << 1.,2.,3.;
  V_t x; x << 2.,3.,4.;
  V_t z; z << 3.,4.,5.;

  pressio::ops::elementwise_multiply(1., x, z, 1., y);
  EXPECT_DOUBLE_EQ( y(0), 7.0);
  EXPECT_DOUBLE_EQ( y(1), 14.0);
  EXPECT_DOUBLE_EQ( y(2), 23.0);

  // test beta=0 with simulated NaN in uninitialized y
  y(0) = std::nan("0");
  pressio::ops::elementwise_multiply(1., x, z, 0., y);
  EXPECT_DOUBLE_EQ( y(0), 6.0);
  EXPECT_DOUBLE_EQ( y(1), 12.0);
  EXPECT_DOUBLE_EQ( y(2), 20.0);
}
