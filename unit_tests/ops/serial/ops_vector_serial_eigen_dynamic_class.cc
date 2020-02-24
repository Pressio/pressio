
#include <gtest/gtest.h>
#include "pressio_ops.hpp"

using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using myvec_t = pressio::containers::Vector<eigvec_t>;

TEST(ops_vector_serial_eigen_dynamic_class,
     resize){

  myvec_t m_v1(11);
  ASSERT_FALSE( m_v1.empty() );
  ASSERT_TRUE( m_v1.extent(0) == 11 );
  pressio::ops::resize(m_v1, 22);
  ASSERT_FALSE( m_v1.empty() );
  ASSERT_TRUE( m_v1.extent(0) == 22 );
}

TEST(ops_vector_serial_eigen_dynamic_class,
     setToScalar){

  myvec_t a(4);
  pressio::ops::fill(a, 1.12);
  EXPECT_DOUBLE_EQ( a(0), 1.12);
  EXPECT_DOUBLE_EQ( a(1), 1.12);
  EXPECT_DOUBLE_EQ( a(2), 1.12);
  EXPECT_DOUBLE_EQ( a(3), 1.12);
}

TEST(ops_vector_serial_eigen_dynamic_class,
     settingZero){

  myvec_t a(4);
  ::pressio::ops::set_zero(a);
  EXPECT_DOUBLE_EQ( a(0), 0.0);
  EXPECT_DOUBLE_EQ( a(1), 0.0);
  EXPECT_DOUBLE_EQ( a(2), 0.0);
  EXPECT_DOUBLE_EQ( a(3), 0.0);
}

TEST(ops_vector_serial_eigen_dynamic_class,
     norm1){

  myvec_t a(4);
  a(0) = 1.; a(1) = 2.;
  a(2) = -1.; a(3) = 3.;
  auto res = pressio::ops::norm1(a);
  EXPECT_DOUBLE_EQ( res, 7.0);
}

TEST(ops_vector_serial_eigen_dynamic_class,
     norm2){

  myvec_t a(3);
  a(0) = 1.; a(1) = 2.; a(2) = -2.;
  auto res = pressio::ops::norm2(a);
  EXPECT_DOUBLE_EQ( res, 3.0);
}

TEST(ops_vector_serial_eigen_dynamic_class,
     minValue){

  myvec_t a(3);
  a(0) = 1.; a(1) = 2.; a(2) = -2.;
  auto res = pressio::ops::min(a);
  EXPECT_DOUBLE_EQ( res, -2.0);
}

TEST(ops_vector_serial_eigen_dynamic_class,
     maxValue){

  myvec_t a(3);
  a(0) = 1.; a(1) = 2.; a(2) = -2.;
  auto res = pressio::ops::max(a);
  EXPECT_DOUBLE_EQ( res, 2.0);
}
