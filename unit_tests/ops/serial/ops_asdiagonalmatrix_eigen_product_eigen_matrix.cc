
#include <gtest/gtest.h>
#include "pressio_ops.hpp"

TEST(ops_asdiagonalmatrix_prod_with, productWithEigenMatrix)
{
  using eigdmat_t = Eigen::MatrixXd;
  using eigdvec_t = Eigen::VectorXd;

  // eigen vector
  pressio::containers::Vector<eigdvec_t> v(6);
  for (int i=0; i<6; ++i) v(i)= (double) i;

  // view v as diagonal matrix
  const auto vD = pressio::containers::asDiagonalMatrix(v);

  //construct matrix
  using mat_t = pressio::containers::DenseMatrix<eigdmat_t>;
  mat_t A(6,3);
  A(0,0) = 1.; A(0,1) = 2.; A(0,2) = 3.;
  A(1,0) = 3.; A(1,1) = 2.; A(1,2) = 1.;
  A(2,0) = 0.; A(2,1) = 0.; A(2,2) = 1.;
  A(3,0) = 0.; A(3,1) = 1.; A(3,2) = 0.;
  A(4,0) = 1.; A(4,1) = 0.; A(4,2) = 0.;
  A(5,0) = 0.; A(5,1) = 1.; A(5,2) = 1.;

  mat_t C(6,3);
  pressio::ops::product(pressio::nontranspose(),
			pressio::nontranspose(),
			1., vD, A, 0., C);

  std::cout << *C.data() << std::endl;
  ASSERT_EQ( C.extent(0), 6 );
  ASSERT_EQ( C.extent(1), 3 );

  EXPECT_DOUBLE_EQ(C(0,0),0.);
  EXPECT_DOUBLE_EQ(C(0,1),0.);
  EXPECT_DOUBLE_EQ(C(0,2),0.);

  EXPECT_DOUBLE_EQ(C(1,0),3.);
  EXPECT_DOUBLE_EQ(C(1,1),2.);
  EXPECT_DOUBLE_EQ(C(1,2),1.);

  EXPECT_DOUBLE_EQ(C(2,0),0.);
  EXPECT_DOUBLE_EQ(C(2,1),0.);
  EXPECT_DOUBLE_EQ(C(2,2),2.);

  EXPECT_DOUBLE_EQ(C(3,0),0.);
  EXPECT_DOUBLE_EQ(C(3,1),3.);
  EXPECT_DOUBLE_EQ(C(3,2),0.);

  EXPECT_DOUBLE_EQ(C(4,0),4.);
  EXPECT_DOUBLE_EQ(C(4,1),0.);
  EXPECT_DOUBLE_EQ(C(4,2),0.);

  EXPECT_DOUBLE_EQ(C(5,0),0.);
  EXPECT_DOUBLE_EQ(C(5,1),5.);
  EXPECT_DOUBLE_EQ(C(5,2),5.);
}
