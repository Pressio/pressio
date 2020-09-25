
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

using em_t = Eigen::MatrixXd;
using w_t = pressio::containers::MultiVector<em_t>;

TEST(containers_multivector_sharedmem_eigen_dyn, Constructor1)
{
  w_t b;
  ASSERT_TRUE( b.data()->data() == nullptr );
  ASSERT_TRUE( b.data()->rows() == 0 );
  ASSERT_TRUE( b.data()->cols() == 0 );
}

TEST(containers_multivector_sharedmem_eigen_dyn, Constructor2)
{
  em_t a (15, 4);
  a.setConstant(1);

  w_t b(a);
  for (auto j=0; j<4; ++j){
    ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 15.);
  }
  // change b should not affect a
  b.data()->setConstant(2.);
  for (auto j=0; j<4; ++j){
    ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 30.);
    ASSERT_EQ(a.col(j).lpNorm<1>(), 15.);
  }

  ASSERT_TRUE( a.data() != b.data()->data() );
  ASSERT_TRUE( a.data() != nullptr );
  ASSERT_TRUE( b.data()->data() != nullptr );
}

TEST(containers_multivector_sharedmem_eigen_dyn, Constructor3)
{
  w_t b(15,4);
  ASSERT_TRUE( b.data()->data() != nullptr );
  ASSERT_TRUE( b.data()->rows() == 15 );
  ASSERT_TRUE( b.data()->cols() == 4 );
}

TEST(containers_multivector_sharedmem_eigen_dyn, Constructor4)
{
  em_t a(15,4);
  a.setConstant(1);
  w_t b (std::move(a));
  for (auto j=0; j<4; ++j){
    ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 15.);
  }

  ASSERT_TRUE( a.data() == nullptr );
  ASSERT_TRUE( b.data() != nullptr );
}

TEST(containers_multivector_sharedmem_eigen_dyn, CopyConstructor)
{
  w_t a(15,4);
  a.data()->setConstant(1.);

  w_t b(a);
  for (auto j=0; j<4; ++j){
    ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 15.);
  }

  // change b should not affect a
  b.data()->setConstant(2.);
  for (auto j=0; j<4; ++j){
    ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 30.);
    ASSERT_EQ(a.data()->col(j).lpNorm<1>(), 15.);
  }

  ASSERT_TRUE( a.data()->data() != b.data()->data() );
  ASSERT_TRUE( a.data()->data() != nullptr );
  ASSERT_TRUE( b.data()->data() != nullptr );
}

TEST(containers_multivector_sharedmem_eigen_dyn, MoveConstructor)
{
  w_t a(15,4);
  a.data()->setConstant(1.);
  const auto ptr = a.data()->data();
  ASSERT_TRUE( ptr != nullptr );

  w_t b(std::move(a));
  for (auto j=0; j<4; ++j){
    ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 15.);
  }

  ASSERT_TRUE( a.data()->data() == nullptr );
  ASSERT_TRUE( b.data()->data() == ptr );
}

TEST(containers_multivector_sharedmem_eigen_dyn, MoveAssign)
{
  w_t b(15,4);
  b.data()->setConstant(1.);
  for (auto j=0; j<4; ++j){
    ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 15.);
  }
  {
    w_t a(15,4);
    const auto ptr = a.data()->data();
    a.data()->setConstant(2.);
    ASSERT_TRUE( ptr != nullptr );
    b = std::move(a);
    for (auto j=0; j<4; ++j){
      ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 30.);
    }

    ASSERT_TRUE( b.data()->data() == ptr );

    // waiting for: https://gitlab.com/libeigen/eigen/-/issues/2000
    //ASSERT_TRUE( a.data()->data() == nullptr );
  }
  ASSERT_TRUE( b.data()->data() != nullptr );
}
