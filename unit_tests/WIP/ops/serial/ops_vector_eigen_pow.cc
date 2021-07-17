
#include <gtest/gtest.h>
#include "pressio_ops.hpp"

TEST(eigenOps, vectorPow)
{
  using vec_t = Eigen::VectorXd;
  vec_t x(6);
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

TEST(eigenOps, vectorAbsPowPos)
{
  using vec_t = Eigen::VectorXd;
  vec_t y(6);
  vec_t x(6);
  for (int i=0; i<6; ++i){
    x(i) = (double) i; x(i)*=-1.;
  }

  ::pressio::ops::abs_pow(y, x, 3.);
  Eigen::VectorXd g(6);
  g(0) = 0.;
  g(1) = 1.;
  g(2) = 8.;
  g(3) = 27.;
  g(4) = 64.;
  g(5) = 125.;
  ASSERT_EQ( y.size(), 6 );
  for (int i=0; i<6; ++i)
    EXPECT_DOUBLE_EQ(y(i), g(i));
}

TEST(eigenOps, vectorAbsPowNeg)
{
  using vec_t = Eigen::VectorXd;
  vec_t y(6);
  vec_t x(6);
  for (int i=0; i<6; ++i){
    x(i) = (double) i; x(i)*=-1.;
  }

  ::pressio::ops::abs_pow(y, x, -3., 0.00001);
  std::cout << y << std::endl;

  Eigen::VectorXd g(6);
  g(0) = 1./0.00001; // because we guard against diving by zero
  g(1) = 1.;
  g(2) = 1./8.;
  g(3) = 1./27.;
  g(4) = 1./64.;
  g(5) = 1./125.;
  ASSERT_EQ( y.size(), 6 );
  for (int i=0; i<6; ++i)
    EXPECT_DOUBLE_EQ(y(i), g(i));
}
