
#include <gtest/gtest.h>
#include "pressio_ops.hpp"

TEST(ops_matrix_matrix_product, eigenDenseDense){
  using namespace pressio;

  using nat_t = Eigen::MatrixXd;
  nat_t A(3,4);
  A << 1.,2.,3.,4., 4.,3.,2.,1., 1.,2.,3.,4.;

  using nat2_t = Eigen::MatrixXd;
  nat2_t B(4,2); 
  B << 1.,2.,3.,3.,4.,4.,2.,2.;

  Eigen::Matrix<double,3,2> C;

  constexpr auto beta  = ::pressio::utils::constants<double>::zero();
  constexpr auto alpha = ::pressio::utils::constants<double>::one();
  ::pressio::ops::product(::pressio::nontranspose(), 
    ::pressio::nontranspose(), alpha, A, B, beta, C);
  EXPECT_DOUBLE_EQ( C(0,0), 27.0);
  EXPECT_DOUBLE_EQ( C(1,0), 23.0);

}//end TEST
