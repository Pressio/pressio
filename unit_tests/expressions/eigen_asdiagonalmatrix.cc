
#include <gtest/gtest.h>
#include "pressio_expressions.hpp"

namespace
{
template <typename T>
void test1(T & a)
{
  {
    const auto D = ::pressio::expressions::asDiagonalMatrix(a);
    EXPECT_EQ( D.extent(0), 6 );
    EXPECT_EQ( D.extent(1), 6 );
    EXPECT_DOUBLE_EQ( D(0,0), 1. );
    EXPECT_DOUBLE_EQ( D(1,1), 5. );
    EXPECT_DOUBLE_EQ( D(5,5), 21. );
  }
}

template <typename T>
void test2(T & a)
{
  {
    // change some entries
    auto D = pressio::expressions::asDiagonalMatrix(a);
    EXPECT_EQ( D.extent(0), 6 );

    // before changing it
    EXPECT_DOUBLE_EQ( D(0,0), 1. );
    EXPECT_DOUBLE_EQ( D(1,1), 5. );
    EXPECT_DOUBLE_EQ( D(2,2), 9. );
    // modify
    D(0,0) = 44.;
    // after
    EXPECT_DOUBLE_EQ( D(0,0), 44. );
    EXPECT_DOUBLE_EQ( D(1,1), 5. );
    EXPECT_DOUBLE_EQ( D(2,2), 9. );
  }
}

template <typename T>
void testConst(const T & a){
  auto D = pressio::expressions::asDiagonalMatrix(a);
  EXPECT_EQ( D.extent(0), 6 );
  EXPECT_DOUBLE_EQ( D(0,0), 44. );
  EXPECT_DOUBLE_EQ( D(1,1), 5. );
  EXPECT_DOUBLE_EQ( D(2,2), 9. );

  // D(0,0) = 1.;
  // does not compile since "a" is const, so should not compile!
}
};

TEST(containers_eigen, asDiagonalMatrix)
{
  using v_t = Eigen::VectorXd;
  v_t a(6);
  a(0) = 1.;
  a(1) = 5.;
  a(2) = 9.;
  a(3) = 13.;
  a(4) = 17.;
  a(5) = 21.;

  test1(a);
  test2(a);
  testConst(a);
}
