
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

using eigdmat_t = Eigen::MatrixXd;
using myMV_t = pressio::containers::MultiVector<eigdmat_t>;


TEST(containers_multi_vector_serial_eigen_dynamic_class,
     productWithEigenVectorViewFromMV){

  //construct multivector
  myMV_t A(6,3);
  A(0,0) = 1.; A(0,1) = 2.; A(0,2) = 3.;
  A(1,0) = 3.; A(1,1) = 2.; A(1,2) = 1.;
  A(2,0) = 0.; A(2,1) = 0.; A(2,2) = 1.;
  A(3,0) = 0.; A(3,1) = 1.; A(3,2) = 0.;
  A(4,0) = 1.; A(4,1) = 0.; A(4,2) = 0.;
  A(5,0) = 0.; A(5,1) = 1.; A(5,2) = 1.;

  // eigen mv
  pressio::containers::MultiVector<Eigen::MatrixXd> b(3,2);
  b(0,0) = 1.; b(0,1) = 2.;
  b(1,0) = 1.; b(1,1) = 2.;
  b(2,0) = 1.; b(2,1) = 2.;

  {
    const auto colVec = pressio::containers::viewColumnVector(b, 0);
    auto c1 = pressio::containers::ops::product(A, colVec);
    ASSERT_EQ( c1.extent(0), 6 );
    EXPECT_DOUBLE_EQ( c1(0), 6.);
    EXPECT_DOUBLE_EQ( c1(1), 6.);
    EXPECT_DOUBLE_EQ( c1(2), 1.);
    EXPECT_DOUBLE_EQ( c1(3), 1.);
    EXPECT_DOUBLE_EQ( c1(4), 1.);
    EXPECT_DOUBLE_EQ( c1(5), 2.);

    pressio::containers::Vector<Eigen::VectorXd> c(6);
    pressio::containers::ops::product(A,colVec,c);
    for (auto i=0; i<c1.extent(0); i++)
      EXPECT_DOUBLE_EQ( c(i), c1(i));
  }

  {
    const auto colVec = pressio::containers::viewColumnVector(b, 1);
    auto c1 = pressio::containers::ops::product(A, colVec);
    ASSERT_EQ( c1.extent(0), 6 );
    EXPECT_DOUBLE_EQ( c1(0), 12.);
    EXPECT_DOUBLE_EQ( c1(1), 12.);
    EXPECT_DOUBLE_EQ( c1(2), 2.);
    EXPECT_DOUBLE_EQ( c1(3), 2.);
    EXPECT_DOUBLE_EQ( c1(4), 2.);
    EXPECT_DOUBLE_EQ( c1(5), 4.);

    pressio::containers::Vector<Eigen::VectorXd> c(6);
    pressio::containers::ops::product(A,colVec,c);
    for (auto i=0; i<c1.extent(0); i++)
      EXPECT_DOUBLE_EQ( c(i), c1(i));
  }
}
