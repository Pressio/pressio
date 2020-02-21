
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using myvec_t = pressio::containers::Vector<eigvec_t>;


TEST(containers_vector_serial_eigen_dynamic_class,
     constructor){

  using vecTrait = pressio::containers::details::traits<myvec_t>;
  ASSERT_TRUE(vecTrait::wrapped_vector_identifier
  == pressio::containers::details::WrappedVectorIdentifier::EigenColDynamic);

  //construct by passing the size
  myvec_t m_v2(5);

  // pass native eigen vector
  eigvec_t e_v1(45);
  e_v1(2) = 2.2; e_v1(4) = 4.4;
  myvec_t m_v1(e_v1);
}

TEST(containers_vector_serial_eigen_dynamic_class,
     constructorAndCheckVals){

  //construct by passing the size
  myvec_t m_v2(5);
  ASSERT_TRUE( m_v2.extent(0) == 5 );
  for (size_t i=0; i<5; i++)
    EXPECT_DOUBLE_EQ( m_v2[i], 0.);

  // pass native eigen vector
  eigvec_t e_v1(45);
  e_v1(2) = 2.2; e_v1(4) = 4.4;
  myvec_t m_v1(e_v1);
  ASSERT_TRUE( m_v1.extent(0) == 45 );
  EXPECT_DOUBLE_EQ( m_v1[2], 2.2);
  EXPECT_DOUBLE_EQ( m_v1[4], 4.4);
}


TEST(containers_vector_serial_eigen_dynamic_class,
     copyConstructor){

  //construct by passing the size
  myvec_t a(3);
  a[0]=1.1; a[1]=1.2; a[2]= 1.3;

  myvec_t b(a);
  EXPECT_DOUBLE_EQ( b[0], 1.1);
  EXPECT_DOUBLE_EQ( b[1], 1.2);
  EXPECT_DOUBLE_EQ( b[2], 1.3);
}


TEST(containers_vector_serial_eigen_dynamic_class,
     assignOp){

  //construct by passing the size
  myvec_t a(3);
  a[0]=1.1; a[1]=1.2; a[2]= 1.3;

  myvec_t b(3);
  b = a;
  EXPECT_DOUBLE_EQ( b[0], 1.1);
  EXPECT_DOUBLE_EQ( b[1], 1.2);
  EXPECT_DOUBLE_EQ( b[2], 1.3);
}


TEST(containers_vector_serial_eigen_dynamic_class,
     queryWrappedData){

  myvec_t m_v1(4);
  ::testing::StaticAssertTypeEq<decltype(m_v1.data()),
				eigvec_t * >();
  const myvec_t m_v2(4);
  ::testing::StaticAssertTypeEq< decltype(m_v2.data()),
				 const eigvec_t * >();
}


TEST(containers_vector_serial_eigen_dynamic_class,
     size){

  myvec_t m_v1(11);
  ASSERT_TRUE( m_v1.extent(0) == 11 );
}

TEST(containers_vector_serial_eigen_dynamic_class,
     resize){

  myvec_t m_v1(11);
  ASSERT_FALSE( m_v1.empty() );
  ASSERT_TRUE( m_v1.extent(0) == 11 );
  pressio::containers::ops::resize(m_v1, 22);
  ASSERT_FALSE( m_v1.empty() );
  ASSERT_TRUE( m_v1.extent(0) == 22 );
}


TEST(containers_vector_serial_eigen_dynamic_class,
     subscriptOperatorSquareBrack){

  myvec_t m_v3(4);
  ::testing::StaticAssertTypeEq<
    decltype(m_v3[1]), double & >();

  ASSERT_TRUE( m_v3.extent(0) == 4 );
  m_v3[0] = 34.0;
  m_v3[1] = 22.5;
  m_v3[2] = 11.5;
  m_v3[3] = 75.0;
  EXPECT_DOUBLE_EQ( m_v3[0], 34.0);
  EXPECT_DOUBLE_EQ( m_v3[1], 22.5);
  EXPECT_DOUBLE_EQ( m_v3[2], 11.5);
  EXPECT_DOUBLE_EQ( m_v3[3], 75.0);
  m_v3[0] = 56.;
  m_v3[3] = 44.;
  EXPECT_DOUBLE_EQ( m_v3[0], 56.0);
  EXPECT_DOUBLE_EQ( m_v3[3], 44.0);

  const myvec_t m_v4(4);
  ::testing::StaticAssertTypeEq<
    decltype(m_v4[1]), const double & >();
}

TEST(containers_vector_serial_eigen_dynamic_class,
     subscriptOperatorParenthesis){

  myvec_t m_v3(4);
  ASSERT_TRUE( m_v3.extent(0) == 4 );
  ::testing::StaticAssertTypeEq<
    decltype(m_v3(1)), double & >();
  m_v3(0) = 34.0;
  m_v3(1) = 22.5;
  m_v3(2) = 11.5;
  m_v3(3) = 75.0;
  EXPECT_DOUBLE_EQ( m_v3(0), 34.0);
  EXPECT_DOUBLE_EQ( m_v3(1), 22.5);
  EXPECT_DOUBLE_EQ( m_v3(2), 11.5);
  EXPECT_DOUBLE_EQ( m_v3(3), 75.0);

  m_v3(0) = 56.;
  m_v3(3) = 44.;
  EXPECT_DOUBLE_EQ( m_v3(0), 56.0);
  EXPECT_DOUBLE_EQ( m_v3(3), 44.0);

  const myvec_t m_v4(4);
  ::testing::StaticAssertTypeEq<
    decltype(m_v4(1)), const double & >();
}


TEST(containers_vector_serial_eigen_dynamic_class,
     matchingSize){

  myvec_t a(4);
  ASSERT_TRUE( a.extent(0) == 4 );
  myvec_t b(6);
  a = b;
  ASSERT_TRUE( a.extent(0) == 6 );
  ASSERT_FALSE( a.extent(0) == 4 );
}


TEST(containers_vector_serial_eigen_dynamic_class,
     setToScalar){

  myvec_t a(4);
  pressio::containers::ops::fill(a, 1.12);
  EXPECT_DOUBLE_EQ( a(0), 1.12);
  EXPECT_DOUBLE_EQ( a(1), 1.12);
  EXPECT_DOUBLE_EQ( a(2), 1.12);
  EXPECT_DOUBLE_EQ( a(3), 1.12);
}


TEST(containers_vector_serial_eigen_dynamic_class,
     assignScalar){

  myvec_t a(4);
  a = 1.12;
  EXPECT_DOUBLE_EQ( a(0), 1.12);
  EXPECT_DOUBLE_EQ( a(1), 1.12);
  EXPECT_DOUBLE_EQ( a(2), 1.12);
  EXPECT_DOUBLE_EQ( a(3), 1.12);
}


TEST(containers_vector_serial_eigen_dynamic_class,
     settingZero){

  myvec_t a(4);
  ::pressio::containers::ops::set_zero(a);
  EXPECT_DOUBLE_EQ( a(0), 0.0);
  EXPECT_DOUBLE_EQ( a(1), 0.0);
  EXPECT_DOUBLE_EQ( a(2), 0.0);
  EXPECT_DOUBLE_EQ( a(3), 0.0);
}

TEST(containers_vector_serial_eigen_dynamic_class,
     norm1){

  myvec_t a(4);
  a(0) = 1.; a(1) = 2.;
  a(2) = -1.; a(3) = 3.;
  auto res = pressio::containers::ops::norm1(a);
  EXPECT_DOUBLE_EQ( res, 7.0);
}

TEST(containers_vector_serial_eigen_dynamic_class,
     norm2){

  myvec_t a(3);
  a(0) = 1.; a(1) = 2.; a(2) = -2.;
  auto res = pressio::containers::ops::norm2(a);
  EXPECT_DOUBLE_EQ( res, 3.0);
}

TEST(containers_vector_serial_eigen_dynamic_class,
     minValue){

  myvec_t a(3);
  a(0) = 1.; a(1) = 2.; a(2) = -2.;
  auto res = pressio::containers::ops::min(a);
  EXPECT_DOUBLE_EQ( res, -2.0);
}

TEST(containers_vector_serial_eigen_dynamic_class,
     maxValue){

  myvec_t a(3);
  a(0) = 1.; a(1) = 2.; a(2) = -2.;
  auto res = pressio::containers::ops::max(a);
  EXPECT_DOUBLE_EQ( res, 2.0);
}


TEST(containers_vector_serial_eigen_dynamic_class,
     CompoundAssignAddOperator){

  myvec_t m_v1(4);
  m_v1[0] = 3.; m_v1[1] = 2.;
  m_v1[2] = 4.; m_v1[3] = 5.;
  myvec_t m_v2(4);
  m_v2[0] = 1.; m_v2[1] = 1.;
  m_v2[2] = 1.; m_v2[3] = 1.;

  m_v1 += m_v2;
  EXPECT_DOUBLE_EQ(m_v1[0], 4.);
  EXPECT_DOUBLE_EQ(m_v1[1], 3.);
  EXPECT_DOUBLE_EQ(m_v1[2], 5.);
  EXPECT_DOUBLE_EQ(m_v1[3], 6.);
}


TEST(containers_vector_serial_eigen_dynamic_class,
     CompoundAssignSubtractOperator){

  myvec_t m_v1(4);
  m_v1[0] = 3.; m_v1[1] = 2.;
  m_v1[2] = 4.; m_v1[3] = 5.;
  myvec_t m_v2(4);
  m_v2[0] = 1.; m_v2[1] = 1.;
  m_v2[2] = 1.; m_v2[3] = 1.;

  m_v1 -= m_v2;
  EXPECT_DOUBLE_EQ(m_v1[0], 2.);
  EXPECT_DOUBLE_EQ(m_v1[1], 1.);
  EXPECT_DOUBLE_EQ(m_v1[2], 3.);
  EXPECT_DOUBLE_EQ(m_v1[3], 4.);
}
