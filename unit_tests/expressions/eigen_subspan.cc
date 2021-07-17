
#include <gtest/gtest.h>
#include "pressio_expressions.hpp"

namespace
{

template <typename T>
void test1(T & A)
{
  {
    const auto sspan = pressio::subspan(A,
						    std::make_pair(0,2),
						    std::make_pair(1,2) );
    EXPECT_EQ( sspan.extent(0), 2 ); EXPECT_EQ( sspan.extent(1), 1 );
    EXPECT_DOUBLE_EQ( sspan(0,0), 2. );
    EXPECT_DOUBLE_EQ( sspan(1,0), 6. );
  }
  {
    const auto sspan = pressio::subspan(A,
						    std::make_pair(0,3),
						    std::make_pair(1,3) );
    EXPECT_EQ( sspan.extent(0), 3 ); EXPECT_EQ( sspan.extent(1), 2 );
    EXPECT_DOUBLE_EQ( sspan(0,0), 2. );  EXPECT_DOUBLE_EQ( sspan(0,1), 3. );
    EXPECT_DOUBLE_EQ( sspan(1,0), 6. );  EXPECT_DOUBLE_EQ( sspan(1,1), 7. );
    EXPECT_DOUBLE_EQ( sspan(2,0), 10. ); EXPECT_DOUBLE_EQ( sspan(2,1), 11. );
  }

  {
    const auto sspan = pressio::subspan(A,
						    std::make_pair(2,4),
						    std::make_pair(1,3) );
    EXPECT_EQ( sspan.extent(0), 2 ); EXPECT_EQ( sspan.extent(1), 2 );
    EXPECT_DOUBLE_EQ( sspan(0,0), 10. );  EXPECT_DOUBLE_EQ( sspan(0,1), 11. );
    EXPECT_DOUBLE_EQ( sspan(1,0), 14. );  EXPECT_DOUBLE_EQ( sspan(1,1), 15. );
  }
  {
    const auto sspan = pressio::subspan(A,
						    std::make_pair(2,3),
						    std::make_pair(0,1) );
    EXPECT_EQ( sspan.extent(0), 1 ); EXPECT_EQ( sspan.extent(1), 1 );
    EXPECT_DOUBLE_EQ( sspan(0,0), 9. );
  }
}

template <typename T>
void test2(T & A)
{
  {
    // change some entries
    auto sspan = pressio::subspan(A,
					      std::make_pair(2,4),
					      std::make_pair(1,3) );
    EXPECT_EQ( sspan.extent(0), 2 ); EXPECT_EQ( sspan.extent(1), 2 );

    // before changing it
    EXPECT_DOUBLE_EQ( sspan(0,0), 10. );  EXPECT_DOUBLE_EQ( sspan(0,1), 11. );
    EXPECT_DOUBLE_EQ( sspan(1,0), 14. );  EXPECT_DOUBLE_EQ( sspan(1,1), 15. );
    // modify
    sspan(0,0) = 44.; sspan(0,1) = 33.;
    // after
    EXPECT_DOUBLE_EQ( sspan(0,0), 44. );  EXPECT_DOUBLE_EQ( sspan(0,1), 33. );
    EXPECT_DOUBLE_EQ( sspan(1,0), 14. );  EXPECT_DOUBLE_EQ( sspan(1,1), 15. );
  }

  {
    // get the native expression
    const auto sspan = pressio::subspan(A,
						    std::make_pair(2,4),
						    std::make_pair(1,3) );
    auto &natEx = *sspan.data();
    EXPECT_EQ( natEx.rows(), 2 ); EXPECT_EQ( natEx.cols(), 2 );
    EXPECT_DOUBLE_EQ( natEx(0,0), 44. );  EXPECT_DOUBLE_EQ( natEx(0,1), 33. );
    EXPECT_DOUBLE_EQ( natEx(1,0), 14. );  EXPECT_DOUBLE_EQ( natEx(1,1), 15. );
  }
}

template <typename T>
void testConst(const T & A)
{
  const auto sspan = pressio::subspan(A,
						  std::make_pair(2,4),
						  std::make_pair(1,3) );
  EXPECT_EQ( sspan.extent(0), 2 ); EXPECT_EQ( sspan.extent(1), 2 );
  EXPECT_DOUBLE_EQ( sspan(0,0), 44. );  EXPECT_DOUBLE_EQ( sspan(0,1), 33. );
  EXPECT_DOUBLE_EQ( sspan(1,0), 14. );  EXPECT_DOUBLE_EQ( sspan(1,1), 15. );

  auto & natEx = *sspan.data();
  EXPECT_DOUBLE_EQ( natEx(0,0), 44. );  EXPECT_DOUBLE_EQ( natEx(0,1), 33. );
  EXPECT_DOUBLE_EQ( natEx(1,0), 14. );  EXPECT_DOUBLE_EQ( natEx(1,1), 15. );
}

};


TEST(matrix_serial_eigen, subspan)
{
  // col-major matrix (which is default in Eigen)
  using eigmat_t = Eigen::MatrixXd;

  eigmat_t A(6,4);
  A(0,0) = 1.;  A(0,1) = 2.;  A(0,2) = 3.;  A(0,3) = 4.;
  A(1,0) = 5.;  A(1,1) = 6.;  A(1,2) = 7.;  A(1,3) = 8.;
  A(2,0) = 9.;  A(2,1) = 10.; A(2,2) = 11.; A(2,3) = 12.;
  A(3,0) = 13.; A(3,1) = 14.; A(3,2) = 15.; A(3,3) = 16.;
  A(4,0) = 17.; A(4,1) = 18.; A(4,2) = 19.; A(4,3) = 20.;
  A(5,0) = 21.; A(5,1) = 22.; A(5,2) = 23.; A(5,3) = 24.;

  test1(A);
  test2(A);
  testConst(A);
}


TEST(matrix_serial_eigen, subspanRowMajor)
{
  using eigmat_t = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;

  eigmat_t A(6,4);
  A(0,0) = 1.;  A(0,1) = 2.;  A(0,2) = 3.;  A(0,3) = 4.;
  A(1,0) = 5.;  A(1,1) = 6.;  A(1,2) = 7.;  A(1,3) = 8.;
  A(2,0) = 9.;  A(2,1) = 10.; A(2,2) = 11.; A(2,3) = 12.;
  A(3,0) = 13.; A(3,1) = 14.; A(3,2) = 15.; A(3,3) = 16.;
  A(4,0) = 17.; A(4,1) = 18.; A(4,2) = 19.; A(4,3) = 20.;
  A(5,0) = 21.; A(5,1) = 22.; A(5,2) = 23.; A(5,3) = 24.;

  test1(A);
  test2(A);
  testConst(A);
}
