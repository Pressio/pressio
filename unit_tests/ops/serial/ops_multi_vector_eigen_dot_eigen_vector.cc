
#include <gtest/gtest.h>
#include "pressio_ops.hpp"

TEST(ops_multi_vector_serial_eigen_dynamic_class,
     dot_WithEigenVector_dynamic){

  using eigdmat_t = Eigen::MatrixXd;
  using myMV_t = pressio::containers::MultiVector<eigdmat_t>;

  //construct multivector
  myMV_t A(6,3);
  A(0,0) = 1.; A(0,1) = 2.; A(0,2) = 3.;
  A(1,0) = 3.; A(1,1) = 2.; A(1,2) = 1.;
  A(2,0) = 0.; A(2,1) = 0.; A(2,2) = 1.;
  A(3,0) = 0.; A(3,1) = 1.; A(3,2) = 0.;
  A(4,0) = 1.; A(4,1) = 0.; A(4,2) = 0.;
  A(5,0) = 0.; A(5,1) = 1.; A(5,2) = 1.;

  // eigen vector
  pressio::containers::Vector<Eigen::VectorXd> b(6);
  b(0) = 1.; b(1) = 1.; b(2) = 1.;
  b(3) = 1.; b(4) = 2.; b(5) = 1.;

  pressio::containers::Vector<Eigen::VectorXd> c(3);
  constexpr auto beta  = ::pressio::utils::constants<double>::zero();
  constexpr auto alpha = ::pressio::utils::constants<double>::one();
  ::pressio::ops::product(::pressio::transpose(), alpha, A, b, beta, c);

  EXPECT_DOUBLE_EQ( c(0), 6.);
  EXPECT_DOUBLE_EQ( c(1), 6.);
  EXPECT_DOUBLE_EQ( c(2), 6.);
}

TEST(ops_multi_vector_serial_eigen_dynamic_class,
     dot_WithEigenVector_static){

  using eigdmat_t = Eigen::MatrixXd;
  using myMV_t = pressio::containers::MultiVector<eigdmat_t>;

  //construct multivector
  myMV_t A(6,3);
  A(0,0) = 1.; A(0,1) = 2.; A(0,2) = 3.;
  A(1,0) = 3.; A(1,1) = 2.; A(1,2) = 1.;
  A(2,0) = 0.; A(2,1) = 0.; A(2,2) = 1.;
  A(3,0) = 0.; A(3,1) = 1.; A(3,2) = 0.;
  A(4,0) = 1.; A(4,1) = 0.; A(4,2) = 0.;
  A(5,0) = 0.; A(5,1) = 1.; A(5,2) = 1.;

  // eigen vector
  pressio::containers::Vector<Eigen::VectorXd> b(6);
  b(0) = 1.; b(1) = 1.; b(2) = 1.;
  b(3) = 1.; b(4) = 2.; b(5) = 1.;

  using eig_v_st = Eigen::Matrix<double, 3, 1>;
  pressio::containers::Vector<eig_v_st> c;
  c.data()->setConstant(0);
  constexpr auto beta  = ::pressio::utils::constants<double>::zero();
  constexpr auto alpha = ::pressio::utils::constants<double>::one();
  ::pressio::ops::product(::pressio::transpose(), alpha, A, b, beta, c);

  ASSERT_EQ(c.extent(0), 3);
  EXPECT_DOUBLE_EQ( c(0), 6.);
  EXPECT_DOUBLE_EQ( c(1), 6.);
  EXPECT_DOUBLE_EQ( c(2), 6.);
}
