
#include <gtest/gtest.h>
#include "pressio_ops.hpp"

TEST(ops_matrix_matrix_product, eigenDenseDense){
  using namespace pressio;

  using nat_t = Eigen::MatrixXd;
  using myA_t = containers::Matrix<nat_t>;
  nat_t a(3,4);
  a << 1.,2.,3.,4., 4.,3.,2.,1., 1.,2.,3.,4.;
  myA_t A(a);
  //  std::cout << *A.data() << "\n";

  using nat2_t = Eigen::MatrixXd;
  using myB_t = containers::Matrix<nat2_t>;
  nat2_t b(4,2); b << 1.,2.,3.,3.,4.,4.,2.,2.;
  myB_t B(b);
  //  std::cout << *B.data();

  static_assert(
   containers::meta::is_dense_matrix_wrapper_eigen<myA_t>::value,"");
  

  containers::Matrix<Eigen::Matrix<double,3,2>> C2;

  constexpr auto beta  = ::pressio::utils::constants::zero<double>();
  constexpr auto alpha = ::pressio::utils::constants::one<double>();
  ::pressio::ops::product(::pressio::nontranspose(), 
    ::pressio::nontranspose(), alpha, A, B, beta, C2);
  EXPECT_DOUBLE_EQ( C2(0,0), 27.0);
  EXPECT_DOUBLE_EQ( C2(1,0), 23.0);

}//end TEST
