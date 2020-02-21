
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

using eigdmat_t = Eigen::MatrixXd;
using myMV_t = pressio::containers::MultiVector<eigdmat_t>;


TEST(containers_multi_vector_serial_eigen_dynamic_class,
     productWithEigenVector){

  //construct multivector
  myMV_t A(6,3);
  A(0,0) = 1.; A(0,1) = 2.; A(0,2) = 3.;
  A(1,0) = 3.; A(1,1) = 2.; A(1,2) = 1.;
  A(2,0) = 0.; A(2,1) = 0.; A(2,2) = 1.;
  A(3,0) = 0.; A(3,1) = 1.; A(3,2) = 0.;
  A(4,0) = 1.; A(4,1) = 0.; A(4,2) = 0.;
  A(5,0) = 0.; A(5,1) = 1.; A(5,2) = 1.;

  // eigen vector
  pressio::containers::Vector<Eigen::VectorXd> b(3);
  b(0) = 1.; b(1) = 1.; b(2) = 1.;

  pressio::containers::Vector<Eigen::VectorXd> c1(6);
  pressio::containers::ops::product(::pressio::nontranspose(), 1., A, b, 0., c1);
  ASSERT_EQ( c1.extent(0), 6 );
  EXPECT_DOUBLE_EQ( c1(0), 6.);
  EXPECT_DOUBLE_EQ( c1(1), 6.);
  EXPECT_DOUBLE_EQ( c1(2), 1.);
  EXPECT_DOUBLE_EQ( c1(3), 1.);
  EXPECT_DOUBLE_EQ( c1(4), 1.);
  EXPECT_DOUBLE_EQ( c1(5), 2.);
}
