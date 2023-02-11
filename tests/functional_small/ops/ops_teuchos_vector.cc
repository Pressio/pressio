
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
