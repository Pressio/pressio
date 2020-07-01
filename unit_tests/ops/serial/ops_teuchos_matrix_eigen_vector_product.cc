
#include <gtest/gtest.h>
#include "pressio_ops.hpp"

TEST(ops_matrix_vector_product, teuchosDenseMatTransposeEigenVec){
  using namespace pressio;

  using natV_t = Eigen::Matrix<double,3,1>;
  natV_t a; a << 1.,2.,3.;
  containers::Vector<natV_t> x(a);

  Teuchos::SerialDenseMatrix<int, double> M(3,2);
  M(0,0) = 1.0; M(0,1) = 2.0;
  M(1,0) = 3.0; M(1,1) = 4.0;
  M(2,0) = 5.0; M(2,1) = 6.0;

  containers::Vector<Eigen::VectorXd> y(2);
  constexpr auto beta  = ::pressio::utils::constants<double>::zero();
  constexpr auto alpha = ::pressio::utils::constants<double>::one();
  ::pressio::ops::product(::pressio::transpose(), alpha, M, x, beta, y);

  EXPECT_DOUBLE_EQ( y[0], 22.0);
  EXPECT_DOUBLE_EQ( y[1], 28.0);
}

TEST(ops_matrix_vector_product, teuchosDenseMatEigenVec){
  using namespace pressio;

  using natV_t = Eigen::Matrix<double,2,1>;
  natV_t a; a << 1.,2.;
  containers::Vector<natV_t> x(a);

  Teuchos::SerialDenseMatrix<int, double> M(3,2);
  M(0,0) = 1.0; M(0,1) = 2.0;
  M(1,0) = 3.0; M(1,1) = 4.0;
  M(2,0) = 5.0; M(2,1) = 6.0;

  containers::Vector<Eigen::VectorXd> y(3);
  constexpr auto beta  = ::pressio::utils::constants<double>::zero();
  constexpr auto alpha = ::pressio::utils::constants<double>::one();
  ::pressio::ops::product(::pressio::nontranspose(), alpha, M, x, beta, y);

  EXPECT_DOUBLE_EQ( y[0], 5.0);
  EXPECT_DOUBLE_EQ( y[1], 11.0);
  EXPECT_DOUBLE_EQ( y[2], 17.0);
}
