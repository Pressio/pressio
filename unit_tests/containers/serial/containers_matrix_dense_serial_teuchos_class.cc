
#include <gtest/gtest.h>
#include "CONTAINERS_MATRIX"

using nat_t = Teuchos::SerialDenseMatrix<int, double>;
using mymat_t = pressio::containers::Matrix<nat_t>;


TEST(containers_matrix_dense_teuchos_class,
     constructor)
{
  // STATIC_ASSERT_IS_CONTAINERS_MATRIX_WRAPPER(mymat_t);
  //  using matTrait = pressio::containers::details::traits<mymat_t>;

  mymat_t m1;
  EXPECT_EQ( m1.extent(0), 0 );
  EXPECT_EQ( m1.extent(1), 0 );

  mymat_t m2(5, 8);
  EXPECT_EQ( m2.extent(0), 5 );
  EXPECT_EQ( m2.extent(1), 8 );

  nat_t eigMat;
  eigMat.shape(56,101);
  mymat_t m3(eigMat);
  EXPECT_EQ( m3.extent(0), 56 );
  EXPECT_EQ( m3.extent(1), 101 );
}


TEST(containers_matrix_dense_teuchos_class,
     queryWrappedData)
{
  mymat_t m1;
  ::testing::StaticAssertTypeEq<decltype(m1.data()),
				nat_t * >();
  const mymat_t m2(45,64);
  ::testing::StaticAssertTypeEq< decltype(m2.data()),
				 const nat_t * >();
}


// TEST(containers_matrix_dense_teuchos_class,
//      sizeResize)
// {
//   mymat_t m1;
//   EXPECT_EQ( m1.extent(0), 0 );
//   EXPECT_EQ( m1.extent(1), 0 );
//   m1.resize(11,45);
//   EXPECT_NE( m1.extent(0), 0 );
//   EXPECT_NE( m1.extent(1), 0 );
//   EXPECT_EQ( m1.extent(0), 11 );
//   EXPECT_EQ( m1.extent(1), 45 );
// }


TEST(containers_matrix_dense_teuchos_class,
     subscriptOperator)
{
  nat_t em1;
  em1.shape(2,3);
  em1(0,0)=34.; em1(0,1)=22.5; em1(0,2)=11.5;
  em1(1,0)=75.; em1(1,1)=3.; em1(1,2)=6.;

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


TEST(containers_matrix_dense_teuchos_class,
     CompoundAssignAddOperator)
{
  mymat_t m1(2,2);
  m1(0,0)=2.; m1(0,1)=4;
  m1(1,0)=3.; m1(1,1)=6.;

  mymat_t m2(2,2);
  m2(0,0)=1.; m2(0,1)=2;
  m2(1,0)=1.; m2(1,1)=2;

  m1 += m2;
  EXPECT_DOUBLE_EQ(m1(0,0), 3.);
  EXPECT_DOUBLE_EQ(m1(0,1), 6.);
  EXPECT_DOUBLE_EQ(m1(1,0), 4.);
  EXPECT_DOUBLE_EQ(m1(1,1), 8.);
}


TEST(containers_matrix_dense_teuchos_class,
     CompoundAssignSubtractOperator)
{
  mymat_t m1(2,2);
  m1(0,0)=2.; m1(0,1)=4;
  m1(1,0)=3.; m1(1,1)=6.;

  mymat_t m2(2,2);
  m2(0,0)=1.; m2(0,1)=2;
  m2(1,0)=1.; m2(1,1)=2;

  m1 -= m2;
  EXPECT_DOUBLE_EQ(m1(0,0), 1.);
  EXPECT_DOUBLE_EQ(m1(0,1), 2.);
  EXPECT_DOUBLE_EQ(m1(1,0), 2.);
  EXPECT_DOUBLE_EQ(m1(1,1), 4.);
}


// TEST(containers_matrix_dense_teuchos_class,
//      addToDiagonal)
// {
//   mymat_t m1(2,2);
//   m1(0,0)=2.; m1(0,1)=4;
//   m1(1,0)=4.; m1(1,1)=6.;

//   m1.addToDiagonal(1.);

//   EXPECT_DOUBLE_EQ(m1(0,0), 3.);
//   EXPECT_DOUBLE_EQ(m1(0,1), 4.);
//   EXPECT_DOUBLE_EQ(m1(1,0), 4.);
//   EXPECT_DOUBLE_EQ(m1(1,1), 7.);
// }
