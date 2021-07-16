
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

using em_t = Eigen::MatrixXd;
using w_t = pressio::containers::MultiVector<em_t>;

TEST(containers_multi_vector_eigen_dynamic, constructorAndCheckVals)
{
  //construct by passing the sizes 
  w_t A(6,3);
  ASSERT_TRUE( A.extent(0) == 6 );
  ASSERT_TRUE( A.numVectors() == 3 );
  for (size_t i=0; i<6; i++)
    for (size_t j=0; j<3; j++)
      EXPECT_DOUBLE_EQ( A(i,j), 0.);
  
  // pass native eigen vector
  em_t eA(45,12);
  eA(2,2) = 2.2;
  eA(4,11) = 4.4;

  w_t B(eA);
  ASSERT_TRUE( B.extent(0) == 45 );
  ASSERT_TRUE( B.numVectors() == 12 );
  ASSERT_FALSE( B.extent(0) == 4 );
  ASSERT_FALSE( B.numVectors() == 1 );
  EXPECT_DOUBLE_EQ( B(2,2), 2.2);
  EXPECT_DOUBLE_EQ( B(4,11), 4.4);
}


TEST(containers_multivector_eigen_dynamic, Constructor1)
{
  w_t b;
  ASSERT_TRUE( b.data()->data() == nullptr );
  ASSERT_TRUE( b.data()->rows() == 0 );
  ASSERT_TRUE( b.data()->cols() == 0 );
}

TEST(containers_multivector_eigen_dynamic, Constructor2)
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

TEST(containers_multivector_eigen_dynamic, Constructor3)
{
  w_t b(15,4);
  ASSERT_TRUE( b.data()->data() != nullptr );
  ASSERT_TRUE( b.data()->rows() == 15 );
  ASSERT_TRUE( b.data()->cols() == 4 );
}

TEST(containers_multivector_eigen_dynamic, Constructor4)
{
  em_t a(15,4);
  a.setConstant(1);
  const auto ptr = a.data();

  w_t b (std::move(a));
  for (auto j=0; j<4; ++j){
    ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 15.);
  }

  ASSERT_TRUE( b.data()->data() == ptr );
  ASSERT_TRUE( a.data() == nullptr );
  ASSERT_TRUE( b.data() != nullptr );
}

TEST(containers_multivector_eigen_dynamic, CopyConstructor)
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

TEST(containers_multivector_eigen_dynamic, MoveConstructor)
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

TEST(containers_multivector_eigen_dynamic, MoveAssign)
{
  w_t b(15,4);
  b.data()->setConstant(1.);
  auto tmp = b.data()->data();
  for (auto j=0; j<4; ++j){
    ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 15.);
  }
  {
    w_t a(15,4);
    tmp = a.data()->data();
    a.data()->setConstant(2.);
    ASSERT_TRUE( tmp != nullptr );
    b = std::move(a);
    for (auto j=0; j<4; ++j){
      ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 30.);
    }

    ASSERT_TRUE( b.data()->data() == tmp );

    // waiting for: https://gitlab.com/libeigen/eigen/-/issues/2000
    //ASSERT_TRUE( a.data()->data() == nullptr );
  }
  ASSERT_TRUE( b.data()->data() == tmp );
  ASSERT_TRUE( b.data()->data() != nullptr );
}
