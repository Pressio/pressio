
#include <gtest/gtest.h>
#include "pressio/ops.hpp"

using V_t = Teuchos::SerialDenseVector<int, double>;

TEST(ops_teuchos, vector_clone)
{
  const int n = 6;
  V_t a(n);
  for (int i = 0; i < 6; ++i){
    a(i) = i + 1.;
  }

  auto b = pressio::ops::clone(a);
  ASSERT_EQ(::pressio::ops::extent(b, 0), 6);
  for (int i = 0; i < 6; ++i){
    ASSERT_DOUBLE_EQ(b(i), a(i));
  }

  // check if b.data() == a.data()
  a(0) = 55.;
  ASSERT_DOUBLE_EQ(b(0), 1.);
}

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
