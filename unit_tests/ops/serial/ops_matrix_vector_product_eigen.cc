
#include <gtest/gtest.h>
#include "pressio_ops.hpp"

TEST(ops_matrix_vector_product, eigenVectorDenseMatrix){
  using namespace pressio;

  using natV_t = Eigen::Matrix<double,3,1>;
  natV_t a; a << 4.,2.,6;
  containers::Vector<natV_t> myV(a);

  using natM_t = Eigen::Matrix<double,3,3>;
  natM_t M;
  M << 1,0,2,2,1,3,0,0,1;
  containers::DenseMatrix<natM_t> myM(M);

  containers::Vector<natV_t> myR;
  constexpr auto beta  = ::pressio::utils::constants<double>::zero();
  constexpr auto alpha = ::pressio::utils::constants<double>::one();
  ::pressio::ops::product(::pressio::nontranspose(), alpha, myM, myV, beta, myR);
  EXPECT_DOUBLE_EQ( myR(0), 16.0);
  EXPECT_DOUBLE_EQ( myR(1), 28.0);
  EXPECT_DOUBLE_EQ( myR(2), 6.0);
}//end TEST
