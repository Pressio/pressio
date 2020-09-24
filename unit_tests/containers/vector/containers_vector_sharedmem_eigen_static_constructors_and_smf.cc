
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

using vec_t = Eigen::Matrix<double, 15, 1>;
using w_t = pressio::containers::Vector<vec_t>;

TEST(containers_vector_sharedmem_eigen_static, Constructor1)
{
  w_t b;
  ASSERT_TRUE( b.extent(0) == 15);
  ASSERT_TRUE( b.data()->data() != nullptr );
}

TEST(containers_vector_sharedmem_eigen_static, Constructor3)
{
  vec_t a;
  a.setConstant(1);
  w_t b(a);
  ASSERT_EQ(b.data()->lpNorm<1>(), 15.);
  ASSERT_TRUE( a.data() != b.data()->data() );
  ASSERT_TRUE( a.data() != nullptr );
  ASSERT_TRUE( b.data()->data() != nullptr );
}

TEST(containers_vector_sharedmem_eigen_static, Constructor4)
{
  /* move semantics for static vectors do not invalidate
     the moved-from object, they just copy the data, see:
     https://gitlab.com/libeigen/eigen/-/issues/2000
  */

  vec_t a;
  a.setConstant(1);
  w_t b(std::move(a));
  ASSERT_EQ(b.data()->lpNorm<1>(), 15.);
  ASSERT_TRUE( a.data() != b.data()->data() );
  ASSERT_TRUE( a.data() != nullptr );
  ASSERT_TRUE( b.data()->data() != nullptr );
}

TEST(containers_vector_sharedmem_eigen_static, CopyConstructor)
{
  w_t a;
  a.data()->setConstant(1.);
  w_t b(a);
  ASSERT_EQ(b.data()->lpNorm<1>(), 15.);
  ASSERT_TRUE( a.data()->data() != b.data()->data() );
  ASSERT_TRUE( a.data()->data() != nullptr );
  ASSERT_TRUE( b.data()->data() != nullptr );
}

// TEST(containers_vector_sharedmem_eigen_static, CopyAssign)
// {
//   w_t a;
//   a.data()->setConstant(1.);
//   w_t b;
//   b = a;
//   ASSERT_EQ(b.data()->lpNorm<1>(), 15.);
//   ASSERT_TRUE( a.data()->data() != b.data()->data() );
//   ASSERT_TRUE( a.data()->data() != nullptr );
//   ASSERT_TRUE( b.data()->data() != nullptr );
// }

TEST(containers_vector_sharedmem_eigen_static, MoveConstructor)
{
  /* move semantics for static vectors do not invalidate
     the moved-from object, they just copy the data, see:
     https://gitlab.com/libeigen/eigen/-/issues/2000
  */

  w_t a;
  a.data()->setConstant(1.);
  const auto ptr = a.data()->data();
  ASSERT_TRUE( ptr != nullptr );
  w_t b(std::move(a));
  ASSERT_EQ(b.data()->lpNorm<1>(), 15.);
  ASSERT_TRUE( a.data()->data() != nullptr );
  ASSERT_TRUE( b.data()->data() != ptr );
}

TEST(containers_vector_sharedmem_eigen_static, MoveAssign)
{
  /* move semantics for static vectors do not invalidate
     the moved-from object, they just copy the data, see:
     https://gitlab.com/libeigen/eigen/-/issues/2000
  */

  w_t a;
  a.data()->setConstant(1.);
  const auto ptr = a.data()->data();
  ASSERT_TRUE( ptr != nullptr );
  w_t b;
  b = std::move(a);
  ASSERT_EQ(b.data()->lpNorm<1>(), 15.);
  ASSERT_TRUE( a.data()->data() != nullptr );
  ASSERT_TRUE( b.data()->data() != ptr );
}
