
#include <gtest/gtest.h>
#include "pressio_ops.hpp"

using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;

TEST(ops_vector_serial_eigen_dynamic_class,
     resize){

  eigvec_t m_v1(11);
  ASSERT_FALSE( m_v1.size() ==0 );
  ASSERT_TRUE( m_v1.size() == 11 );
  pressio::ops::resize(m_v1, 22);
  ASSERT_FALSE( m_v1.size()==0 );
  ASSERT_TRUE( m_v1.size() == 22 );
}

TEST(ops_vector_serial_eigen_dynamic_class,
     setToScalar){

  eigvec_t a(4);
  pressio::ops::fill(a, 1.12);
  EXPECT_DOUBLE_EQ( a(0), 1.12);
  EXPECT_DOUBLE_EQ( a(1), 1.12);
  EXPECT_DOUBLE_EQ( a(2), 1.12);
  EXPECT_DOUBLE_EQ( a(3), 1.12);
}

TEST(ops_vector_serial_eigen_dynamic_class,
     settingZero){

  eigvec_t a(4);
  ::pressio::ops::set_zero(a);
  EXPECT_DOUBLE_EQ( a(0), 0.0);
  EXPECT_DOUBLE_EQ( a(1), 0.0);
  EXPECT_DOUBLE_EQ( a(2), 0.0);
  EXPECT_DOUBLE_EQ( a(3), 0.0);
}

TEST(ops_vector_serial_eigen_dynamic_class,
     norm1){

  eigvec_t a(4);
  a(0) = 1.; a(1) = 2.;
  a(2) = -1.; a(3) = 3.;
  auto res = pressio::ops::norm1(a);
  EXPECT_DOUBLE_EQ( res, 7.0);
}

TEST(ops_vector_serial_eigen_dynamic_class,
     norm2){

  eigvec_t a(3);
  a(0) = 1.; a(1) = 2.; a(2) = -2.;
  auto res = pressio::ops::norm2(a);
  EXPECT_DOUBLE_EQ( res, 3.0);
}

TEST(ops_vector_serial_eigen_dynamic_class,
     minValue){

  eigvec_t a(3);
  a(0) = 1.; a(1) = 2.; a(2) = -2.;
  auto res = pressio::ops::min(a);
  EXPECT_DOUBLE_EQ( res, -2.0);
}

TEST(ops_vector_serial_eigen_dynamic_class,
     maxValue){

  eigvec_t a(3);
  a(0) = 1.; a(1) = 2.; a(2) = -2.;
  auto res = pressio::ops::max(a);
  EXPECT_DOUBLE_EQ( res, 2.0);
}
