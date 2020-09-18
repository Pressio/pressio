
#include <gtest/gtest.h>
#include "pressio_ops.hpp"

using nat_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using mymat_t = pressio::containers::DenseMatrix<nat_t>;

TEST(ops_matrix_dense_eigen_dynamic_class,
     sizeResize)
{
  mymat_t m1;
  EXPECT_EQ( m1.extent(0), 0 );
  EXPECT_EQ( m1.extent(1), 0 );
  pressio::ops::resize(m1, 11, 45);
  EXPECT_NE( m1.extent(0), 0 );
  EXPECT_NE( m1.extent(1), 0 );
  EXPECT_EQ( m1.extent(0), 11 );
  EXPECT_EQ( m1.extent(1), 45 );
}


// TEST(ops_matrix_dense_eigen_dynamic_class,
//      addToDiagonal)
// {
//   nat_t em1;
//   em1.resize(2,2);
//   em1 << 2., 4., 3., 6.;
//   mymat_t m1(em1);

//   m1.addToDiagonal(1.);

//   EXPECT_DOUBLE_EQ(m1(0,0), 3.);
//   EXPECT_DOUBLE_EQ(m1(0,1), 4.);
//   EXPECT_DOUBLE_EQ(m1(1,0), 3.);
//   EXPECT_DOUBLE_EQ(m1(1,1), 7.);
// }

// TEST(ops_matrix_dense_eigen_dynamic_class, additionOperator)
// {
//   using nat_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
//   using mymat_t = containers::DenseMatrix<nat_t>;

//   nat_t em1;
//   em1.resize(2,2);
//   em1 << 2., 4., 3., 6.;
//   mymat_t m1(em1);

//   nat_t em2;
//   em2.resize(2,2);
//   em2 << 1., 2., 1., 2.;
//   mymat_t m2(em2);

//   mymat_t res = m1 + m2;
//   EXPECT_DOUBLE_EQ(res(0,0), 3.);
//   EXPECT_DOUBLE_EQ(res(0,1), 6.);
//   EXPECT_DOUBLE_EQ(res(1,0), 4.);
//   EXPECT_DOUBLE_EQ(res(1,1), 8.);
// }

// TEST(ops_matrix_dense_eigen_dynamic_class, substractOperator)
// {
//   using nat_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
//   using mymat_t = containers::DenseMatrix<nat_t>;

//   nat_t em1;
//   em1.resize(2,2);
//   em1 << 2., 4., 3., 6.;
//   mymat_t m1(em1);

//   nat_t em2;
//   em2.resize(2,2);
//   em2 << 1., 2., 1., 2.;
//   mymat_t m2(em2);

//   mymat_t res = m1 - m2;
//   EXPECT_DOUBLE_EQ(res(0,0), 1.);
//   EXPECT_DOUBLE_EQ(res(0,1), 2.);
//   EXPECT_DOUBLE_EQ(res(1,0), 2.);
//   EXPECT_DOUBLE_EQ(res(1,1), 4.);
// }


// TEST(ops_matrix_dense_eigen_dynamic_class, starOperator)
// {
//   using nat_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
//   using mymat_t = containers::DenseMatrix<nat_t>;

//   nat_t em1;
//   em1.resize(2,2);
//   em1 << 2., 4., 3., 6.;
//   mymat_t m1(em1);

//   nat_t em2;
//   em2.resize(2,2);
//   em2 << 1., 2., 1., 2.;
//   mymat_t m2(em2);

//   mymat_t res = m1 * m2;
//   EXPECT_DOUBLE_EQ(res(0,0), 6.);
//   EXPECT_DOUBLE_EQ(res(0,1), 12.);
//   EXPECT_DOUBLE_EQ(res(1,0), 9.);
//   EXPECT_DOUBLE_EQ(res(1,1), 18.);
// }



// TEST(ops_matrix_dense_eigen_dynamic_class, transposeDynamic)
// {
//   using nat_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
//   using mymat_t = containers::DenseMatrix<nat_t>;

//   nat_t em1;
//   em1.resize(2,3);
//   em1 << 34.0, 22.5, 11.5, 75., 3., 6.;
//   std::cout << em1 << std::endl;

//   mymat_t m1(em1);
//   EXPECT_DOUBLE_EQ( m1(0,0), 34.0);
//   EXPECT_DOUBLE_EQ( m1(0,1), 22.5);
//   EXPECT_DOUBLE_EQ( m1(0,2), 11.5);
//   EXPECT_DOUBLE_EQ( m1(1,0), 75.);
//   EXPECT_DOUBLE_EQ( m1(1,1), 3.0);
//   EXPECT_DOUBLE_EQ( m1(1,2), 6.0);

//   auto tm1 = containers::mat_ops::transpose(m1);
//   std::cout << *tm1.data() << std::endl;
//   EXPECT_DOUBLE_EQ( tm1(0,0), 34.0);
//   EXPECT_DOUBLE_EQ( tm1(1,0), 22.5);
//   EXPECT_DOUBLE_EQ( tm1(2,0), 11.5);
//   EXPECT_DOUBLE_EQ( tm1(0,1), 75.);
//   EXPECT_DOUBLE_EQ( tm1(1,1), 3.0);
//   EXPECT_DOUBLE_EQ( tm1(2,1), 6.0);
// }


// TEST(ops_matrix_dense_eigen_dynamic_class, transposeStatic)
// {
//   using nat_t = Eigen::Matrix<double, 2, 3>;
//   using mymat_t = containers::DenseMatrix<nat_t>;

//   nat_t em1;
//   em1 << 34.0, 22.5, 11.5, 75., 3., 6.;
//   //std::cout << em1 << std::endl;

//   mymat_t m1(em1);
//   EXPECT_DOUBLE_EQ( m1(0,0), 34.0);
//   EXPECT_DOUBLE_EQ( m1(0,1), 22.5);
//   EXPECT_DOUBLE_EQ( m1(0,2), 11.5);
//   EXPECT_DOUBLE_EQ( m1(1,0), 75.);
//   EXPECT_DOUBLE_EQ( m1(1,1), 3.0);
//   EXPECT_DOUBLE_EQ( m1(1,2), 6.0);

//   auto tm1 = containers::mat_ops::transpose(m1);
//   //std::cout << *tm1.data() << std::endl;
//   EXPECT_DOUBLE_EQ( tm1(0,0), 34.0);
//   EXPECT_DOUBLE_EQ( tm1(1,0), 22.5);
//   EXPECT_DOUBLE_EQ( tm1(2,0), 11.5);
//   EXPECT_DOUBLE_EQ( tm1(0,1), 75.);
//   EXPECT_DOUBLE_EQ( tm1(1,1), 3.0);
//   EXPECT_DOUBLE_EQ( tm1(2,1), 6.0);
// }
