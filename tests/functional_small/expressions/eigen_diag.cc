
#include <gtest/gtest.h>
#include "pressio/expressions.hpp"

namespace{

  template <typename T>
  void test1(T & A)
  {
    {
      const auto diagvals = pressio::diag(A);
      EXPECT_EQ( diagvals.extent(), 4 );
      EXPECT_DOUBLE_EQ( diagvals(0), 1.2 );
      EXPECT_DOUBLE_EQ( diagvals(1), 6.2 );
      EXPECT_DOUBLE_EQ( diagvals(2), 11.2 );
    }
  }

  template <typename T>
  void test2(T & A)
  {
    {
      // change some entries
      auto diagvals = pressio::diag(A);
      EXPECT_EQ( diagvals.extent(), 4 );
      // before changing it
      EXPECT_DOUBLE_EQ( diagvals(0), 1.2 );
      EXPECT_DOUBLE_EQ( diagvals(1), 6.2 );
      EXPECT_DOUBLE_EQ( diagvals(2), 11.2 );
      // modify
      diagvals(0) = 44.;
      diagvals(1) = 6.;
      // after
      EXPECT_DOUBLE_EQ( diagvals(0), 44. );  
      EXPECT_DOUBLE_EQ( diagvals(1), 6. );  
    }

    {
      // get the native expression
      const auto diagvals = pressio::diag(A);
      auto &natEx = diagvals.native();
      EXPECT_DOUBLE_EQ( natEx(0), 44. );  
      EXPECT_DOUBLE_EQ( natEx(1), 6. ); 
    }
  }

  template <typename T>
  void testConst(const T & A){
    const  auto diagvals = pressio::diag(A);
    EXPECT_EQ( diagvals.extent(), 4 );
    EXPECT_DOUBLE_EQ( diagvals(0), 44. );
    EXPECT_DOUBLE_EQ( diagvals(1), 6. );
    EXPECT_DOUBLE_EQ( diagvals(2), 11.2 );

    auto & natEx = diagvals.native();
    EXPECT_DOUBLE_EQ( natEx(0), 44. );  
    EXPECT_DOUBLE_EQ( natEx(1), 6. );
  }
};

TEST(expressions_eigen, diag)
{
  // col-major matrix (which is default in Eigen)
  using eigmat_t = Eigen::MatrixXd;

  eigmat_t A(4,4);
  A(0,0) = 1.2;  A(0,1) = 2.;  A(0,2) = 3.;  A(0,3) = 4.;
  A(1,0) = 5.;  A(1,1) = 6.2;  A(1,2) = 7.;  A(1,3) = 8.;
  A(2,0) = 9.;  A(2,1) = 10.; A(2,2) = 11.2; A(2,3) = 12.;
  A(3,0) = 13.; A(3,1) = 14.; A(3,2) = 15.; A(3,3) = 16.;
  test1(A);
  test2(A);
  testConst(A);
}

TEST(expressions_eigen, diagRowMajor)
{
  using eigmat_t = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;
  eigmat_t A(4,4);
  A(0,0) = 1.2;  A(0,1) = 2.;  A(0,2) = 3.;  A(0,3) = 4.;
  A(1,0) = 5.;  A(1,1) = 6.2;  A(1,2) = 7.;  A(1,3) = 8.;
  A(2,0) = 9.;  A(2,1) = 10.; A(2,2) = 11.2; A(2,3) = 12.;
  A(3,0) = 13.; A(3,1) = 14.; A(3,2) = 15.; A(3,3) = 16.;
  test1(A);
  test2(A);
  testConst(A);
}
