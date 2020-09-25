
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

using vec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using w_t = pressio::containers::Vector<vec_t>;

TEST(containers_vector_sharedmem_eigen_dyn, Constructor1)
{
  w_t b;
  ASSERT_TRUE( b.data()->data() == nullptr );
  ASSERT_TRUE( b.data()->size() == 0 );
}

TEST(containers_vector_sharedmem_eigen_dyn, Constructor2)
{
  vec_t a (15);
  a.setConstant(1);

  w_t b(a);
  ASSERT_EQ(b.data()->lpNorm<1>(), 15.);
  // change b should not affect a
  b.data()->setConstant(2.);
  ASSERT_EQ(b.data()->lpNorm<1>(), 30.);
  ASSERT_EQ(a.lpNorm<1>(), 15.);

  ASSERT_TRUE( a.data() != b.data()->data() );
  ASSERT_TRUE( a.data() != nullptr );
  ASSERT_TRUE( b.data()->data() != nullptr );
}

TEST(containers_vector_sharedmem_eigen_dyn, Constructor3)
{
  w_t b(15);
  ASSERT_TRUE( b.data()->data() != nullptr );
  ASSERT_TRUE( b.data()->size() == 15 );
}

TEST(containers_vector_sharedmem_eigen_dyn, Constructor4)
{
  vec_t a(15);
  a.setConstant(1);
  w_t b (std::move(a));
  ASSERT_EQ(b.data()->lpNorm<1>(), 15.);
  ASSERT_TRUE( a.data() == nullptr );
  ASSERT_TRUE( b.data() != nullptr );
}

TEST(containers_vector_sharedmem_eigen_dyn, CopyConstructor)
{
  w_t a(15);
  a.data()->setConstant(1.);

  w_t b(a);
  ASSERT_EQ(b.data()->lpNorm<1>(), 15.);
  // change b should not affect a
  b.data()->setConstant(2.);
  ASSERT_EQ(b.data()->lpNorm<1>(), 30.);
  ASSERT_EQ(a.data()->lpNorm<1>(), 15.);

  ASSERT_TRUE( a.data()->data() != b.data()->data() );
  ASSERT_TRUE( a.data()->data() != nullptr );
  ASSERT_TRUE( b.data()->data() != nullptr );
}

TEST(containers_vector_sharedmem_eigen_dyn, MoveConstructor)
{
  w_t a(15);
  a.data()->setConstant(1.);
  const auto ptr = a.data()->data();
  ASSERT_TRUE( ptr != nullptr );

  w_t b(std::move(a));
  ASSERT_EQ(b.data()->lpNorm<1>(), 15.);
  ASSERT_TRUE( a.data()->data() == nullptr );
  ASSERT_TRUE( b.data()->data() == ptr );
}

TEST(containers_vector_sharedmem_eigen_dyn, MoveAssign)
{
  w_t b(15);
  b.data()->setConstant(1.);
  ASSERT_EQ(b.data()->lpNorm<1>(), 15.);
  {
    w_t a(15);
    const auto ptr = a.data()->data();
    a.data()->setConstant(2.);
    ASSERT_TRUE( ptr != nullptr );
    b = std::move(a);
    ASSERT_EQ(b.data()->lpNorm<1>(), 30.);

    ASSERT_TRUE( b.data()->data() == ptr );

    // waiting for: https://gitlab.com/libeigen/eigen/-/issues/2000
    //ASSERT_TRUE( ptr == nullptr );
  }
  ASSERT_TRUE( b.data()->data() != nullptr );
}


TEST(containers_vector_sharedmem_eigen_dyn,
     queryWrappedData){

  w_t m_v1(4);
  ::testing::StaticAssertTypeEq<decltype(m_v1.data()),
				vec_t * >();
  const w_t m_v2(4);
  ::testing::StaticAssertTypeEq< decltype(m_v2.data()),
				 const vec_t * >();
}


TEST(containers_vector_sharedmem_eigen_dyn,
     size){

  w_t m_v1(11);
  ASSERT_TRUE( m_v1.extent(0) == 11 );
}


TEST(containers_vector_sharedmem_eigen_dyn,
     subscriptOperatorSquareBrack){

  w_t m_v3(4);
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

  const w_t m_v4(4);
  ::testing::StaticAssertTypeEq<
    decltype(m_v4[1]), const double & >();
}

TEST(containers_vector_sharedmem_eigen_dyn,
     subscriptOperatorParenthesis){

  w_t m_v3(4);
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

  const w_t m_v4(4);
  ::testing::StaticAssertTypeEq<
    decltype(m_v4(1)), const double & >();
}


TEST(containers_vector_sharedmem_eigen_dyn,
     matchingSize){

  w_t a(4);
  ASSERT_TRUE( a.extent(0) == 4 );
  // w_t b(6);
  // a = b;
  // ASSERT_TRUE( a.extent(0) == 6 );
  // ASSERT_FALSE( a.extent(0) == 4 );
}

