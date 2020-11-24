
#include <gtest/gtest.h>
#include "pressio_ops.hpp"

TEST(eigenOps, vectorAbs)
{
  using vec_t = Eigen::VectorXd;

  pressio::containers::Vector<vec_t> x(6);
  for (int i=0; i<6; ++i) x(i) = - (double) i;

  pressio::containers::Vector<vec_t> y(6);
  ::pressio::ops::abs(y,x);

  Eigen::VectorXd g(6);
  g(0) = 0.;
  g(1) = 1.;
  g(2) = 2.;
  g(3) = 3.;
  g(4) = 4.;
  g(5) = 5.;
  ASSERT_EQ( y.extent(0), 6 );
  for (int i=0; i<6; ++i) EXPECT_DOUBLE_EQ(y(i), g(i));
}
