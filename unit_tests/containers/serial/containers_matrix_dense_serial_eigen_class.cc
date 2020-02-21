
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

using nat_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using mymat_t = pressio::containers::Matrix<nat_t>;


TEST(containers_matrix_dense_eigen_dynamic_class,
     constructor)
{
  mymat_t m1;
  EXPECT_EQ( m1.extent(0), 0 );
  EXPECT_EQ( m1.extent(1), 0 );

  mymat_t m2(5, 8);
  EXPECT_EQ( m2.extent(0), 5 );
  EXPECT_EQ( m2.extent(1), 8 );

  nat_t eigMat;
  eigMat.resize(56,101);
  mymat_t m3(eigMat);
  EXPECT_EQ( m3.extent(0), 56 );
  EXPECT_EQ( m3.extent(1), 101 );
}


TEST(containers_matrix_dense_eigen_dynamic_class,
     queryWrappedData)
{
  mymat_t m1;
  ::testing::StaticAssertTypeEq<decltype(m1.data()),
				nat_t * >();
  const mymat_t m2(45,64);
  ::testing::StaticAssertTypeEq< decltype(m2.data()),
				 const nat_t * >();
}

TEST(containers_matrix_dense_eigen_dynamic_class,
     subscriptOperator)
{
  nat_t em1;
  em1.resize(2,3);
  em1 << 34.0, 22.5, 11.5, 75., 3., 6.;

  mymat_t m1(em1);
  EXPECT_DOUBLE_EQ( m1(0,0), 34.0);
  EXPECT_DOUBLE_EQ( m1(0,1), 22.5);
  EXPECT_DOUBLE_EQ( m1(0,2), 11.5);
  EXPECT_DOUBLE_EQ( m1(1,0), 75.);
  EXPECT_DOUBLE_EQ( m1(1,1), 3.0);
  EXPECT_DOUBLE_EQ( m1(1,2), 6.0);

  mymat_t m3(em1);
  EXPECT_EQ( m3.extent(0), 2 );
  EXPECT_EQ( m3.extent(1), 3 );
  m3(1,1) = 55.;
  m3(0,2) = -12.;
  EXPECT_DOUBLE_EQ( m3(1,1), 55.);
  EXPECT_DOUBLE_EQ( m3(0,2), -12.);
}



TEST(containers_matrix_dense_eigen_dynamic_class,
     CompoundAssignAddOperator)
{
  nat_t em1;
  em1.resize(2,2);
  em1 << 2., 4., 3., 6.;
  mymat_t m1(em1);

  nat_t em2;
  em2.resize(2,2);
  em2 << 1., 2., 1., 2.;
  mymat_t m2(em2);

  m1 += m2;
  EXPECT_DOUBLE_EQ(m1(0,0), 3.);
  EXPECT_DOUBLE_EQ(m1(0,1), 6.);
  EXPECT_DOUBLE_EQ(m1(1,0), 4.);
  EXPECT_DOUBLE_EQ(m1(1,1), 8.);
}


TEST(containers_matrix_dense_eigen_dynamic_class,
     CompoundAssignSubtractOperator)
{
  nat_t em1;
  em1.resize(2,2);
  em1 << 2., 4., 3., 6.;
  mymat_t m1(em1);

  nat_t em2;
  em2.resize(2,2);
  em2 << 1., 2., 1., 2.;
  mymat_t m2(em2);

  m1 -= m2;
  EXPECT_DOUBLE_EQ(m1(0,0), 1.);
  EXPECT_DOUBLE_EQ(m1(0,1), 2.);
  EXPECT_DOUBLE_EQ(m1(1,0), 2.);
  EXPECT_DOUBLE_EQ(m1(1,1), 4.);
}

