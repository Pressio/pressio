
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

using vec_t = Eigen::VectorXd;
using w_t = pressio::containers::Vector<vec_t>;

TEST(containers_vector_sharedmem_eigen_dyn, Constructor1)
{
  w_t b;
  ASSERT_TRUE( b.data()->data() == nullptr );
}

TEST(containers_vector_sharedmem_eigen_dyn, Constructor2)
{
  vec_t a (15);
  a.setConstant(1);
  w_t b(a);
  ASSERT_EQ(b.data()->lpNorm<1>(), 15.);
  ASSERT_TRUE( a.data() != b.data()->data() );
  ASSERT_TRUE( a.data() != nullptr );
  ASSERT_TRUE( b.data()->data() != nullptr );
}

TEST(containers_vector_sharedmem_eigen_dyn, Constructor3)
{
  vec_t a(15);
  a.setConstant(1);
  w_t b (std::move(a));
  ASSERT_EQ(b.data()->lpNorm<1>(), 15.);
  ASSERT_TRUE( a.data() == nullptr );
  ASSERT_TRUE( b.data() != nullptr );
}

TEST(containers_vector_sharedmem_eigen_dyn, Constructor4)
{
  w_t b(15);
  ASSERT_TRUE( b.data()->data() != nullptr );
}

TEST(containers_vector_sharedmem_eigen_dyn, CopyConstructor)
{
  w_t a(15);
  a.data()->setConstant(1.);
  w_t b(a);
  ASSERT_EQ(b.data()->lpNorm<1>(), 15.);
  ASSERT_TRUE( a.data()->data() != b.data()->data() );
  ASSERT_TRUE( a.data()->data() != nullptr );
  ASSERT_TRUE( b.data()->data() != nullptr );
}

TEST(containers_vector_sharedmem_eigen_dyn, CopyAssign)
{
  w_t a(15);
  a.data()->setConstant(1.);
  w_t b(15);
  b = a;
  ASSERT_EQ(b.data()->lpNorm<1>(), 15.);
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

// TEST(containers_vector_sharedmem_eigen_dyn, MoveAssign)
// {
//   // waiting for: https://gitlab.com/libeigen/eigen/-/issues/2000

//   // vec_t a(3);
//   // const auto ptr = a.data();
//   // ASSERT_TRUE( ptr != nullptr );
//   // vec_t b(3);
//   // b = std::move(a);
//   // ASSERT_TRUE( a.data() == nullptr );


//   // w_t a(3);
//   // const auto ptr = a.data()->data();
//   // ASSERT_TRUE( ptr != nullptr );
//   // w_t b(3);
//   // b = std::move(a);
//   // ASSERT_TRUE( a.data()->data() == nullptr );
//   // //ASSERT_TRUE( b.data()->data() == ptr );
// }
