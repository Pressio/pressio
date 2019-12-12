
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"

using eigdmat_t = Eigen::MatrixXd;
using myMV_t = pressio::containers::MultiVector<eigdmat_t>;

namespace{

  void test1(myMV_t & A)
  {
    {
      auto c1 = A.viewColumnVector(0);
      c1[0] = 0.5;
      ASSERT_EQ( c1.size(), 6 );
      EXPECT_DOUBLE_EQ( c1[0], .5);
      EXPECT_DOUBLE_EQ( c1[1], 3.);
      EXPECT_DOUBLE_EQ( c1[2], 0.);
      EXPECT_DOUBLE_EQ( c1[3], 0.);
      EXPECT_DOUBLE_EQ( c1[4], 1.);
      EXPECT_DOUBLE_EQ( c1[5], 0.);
    }
    {
      auto c1 = A.viewColumnVector(0);
      auto natExpr = c1();
      natExpr[0] = 0.35;
      EXPECT_DOUBLE_EQ( natExpr[0], .35);
      EXPECT_DOUBLE_EQ( natExpr[1], 3.);
      EXPECT_DOUBLE_EQ( natExpr[2], 0.);
      EXPECT_DOUBLE_EQ( natExpr[3], 0.);
      EXPECT_DOUBLE_EQ( natExpr[4], 1.);
      EXPECT_DOUBLE_EQ( natExpr[5], 0.);
    }
  }

  void test2(const myMV_t & A)
  {
    {
      const auto c1 = A.viewColumnVector(0);
      ASSERT_EQ( c1.size(), 6 );
      EXPECT_DOUBLE_EQ( c1[0], .35);
      EXPECT_DOUBLE_EQ( c1[1], 3.);
      EXPECT_DOUBLE_EQ( c1[2], 0.);
      EXPECT_DOUBLE_EQ( c1[3], 0.);
      EXPECT_DOUBLE_EQ( c1[4], 1.);
      EXPECT_DOUBLE_EQ( c1[5], 0.);
    }
    {
      const auto c1 = A.viewColumnVector(0);
      const auto natExpr = c1();
      EXPECT_DOUBLE_EQ( natExpr[0], .35);
      EXPECT_DOUBLE_EQ( natExpr[1], 3.);
      EXPECT_DOUBLE_EQ( natExpr[2], 0.);
      EXPECT_DOUBLE_EQ( natExpr[3], 0.);
      EXPECT_DOUBLE_EQ( natExpr[4], 1.);
      EXPECT_DOUBLE_EQ( natExpr[5], 0.);
    }
  }
}

TEST(containers_multi_vector_serial_eigen_dynamic_class,
     EigenVectorViewFromMV){

  //construct multivector
  myMV_t A(6,3);
  A(0,0) = 1.; A(0,1) = 2.; A(0,2) = 3.;
  A(1,0) = 3.; A(1,1) = 2.; A(1,2) = 1.;
  A(2,0) = 0.; A(2,1) = 0.; A(2,2) = 1.;
  A(3,0) = 0.; A(3,1) = 1.; A(3,2) = 0.;
  A(4,0) = 1.; A(4,1) = 0.; A(4,2) = 0.;
  A(5,0) = 0.; A(5,1) = 1.; A(5,2) = 1.;

  test1(A);
  test2(A);
}
