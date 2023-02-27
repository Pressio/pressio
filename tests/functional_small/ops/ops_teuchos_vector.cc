
#include <gtest/gtest.h>
#include "pressio/ops.hpp"

using V_t = Teuchos::SerialDenseVector<int, double>;

TEST(ops_teuchos, vector_scale)
{
  const int n = 6;
  V_t a(n);
  for (int i = 0; i < 6; ++i){
    a(i) = i + 1.;
  }

  const double s = 3.;
  pressio::ops::scale(a, s);
  for (int i=0; i<6; ++i){
    ASSERT_DOUBLE_EQ(a(i), (i + 1.) * s);
  }

  // check NaN injection with zero scaling
  a(0) = std::nan("0");
  pressio::ops::scale(a, 0.);
  EXPECT_DOUBLE_EQ(a(0), 0.);
}

TEST(ops_teuchos, vector_norms)
{
  const int n = 2;
  V_t a(n);
  a(0) = -3.;
  a(1) = -4.;
  ASSERT_DOUBLE_EQ(pressio::ops::norm1(a), std::abs(-3) + std::abs(-4));
  ASSERT_DOUBLE_EQ(pressio::ops::norm2(a), std::sqrt( 9. + 16. ));
}
