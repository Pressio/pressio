
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

#ifdef USE_COLMAJ
  using em_t = Eigen::Matrix<double, 6, 4, Eigen::ColMajor>;
#elif USE_ROWMAJ
  using em_t = Eigen::Matrix<double, 6, 4, Eigen::RowMajor>;
#endif

using w_t = pressio::containers::DenseMatrix<em_t>;

TEST(containers_dense_matrix_eigen_static, checkStorage)
{
  //construct by passing the sizes 
  em_t A;
  ASSERT_TRUE( A.rows() == 6 );
  ASSERT_TRUE( A.cols() == 4 );
  for (size_t j=0; j<A.cols(); j++)
    A.col(j).setConstant(j);

#ifdef USE_COLMAJ
  ASSERT_EQ( *A.data(), 0.0);
  ASSERT_EQ( *(A.data()+1), 0.0);
#elif USE_ROWMAJ
  ASSERT_EQ( *A.data(), 0.0);
  ASSERT_EQ( *(A.data()+1), 1.0);
#endif
}

TEST(containers_dense_matrix_eigen_static, Constructor1)
{
  w_t b;
  ASSERT_TRUE( b.data()->data() != nullptr );
  ASSERT_TRUE( b.data()->rows() == 6 );
  ASSERT_TRUE( b.data()->cols() == 4 );
}

TEST(containers_dense_matrix_eigen_static, Constructor2)
{
  em_t a;
  a.setConstant(1);

  w_t b(a);
  for (auto j=0; j<4; ++j){
    ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 6.);
  }
  // change b should not affect a
  b.data()->setConstant(2.);
  for (auto j=0; j<4; ++j){
    ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 12.);
    ASSERT_EQ(a.col(j).lpNorm<1>(), 6.);
  }

  ASSERT_TRUE( a.data() != b.data()->data() );
  ASSERT_TRUE( a.data() != nullptr );
  ASSERT_TRUE( b.data()->data() != nullptr );
}

TEST(containers_dense_matrix_eigen_static, Constructor3)
{
  /* move semantics for static do not invalidate
     the moved-from object, they just copy the data, see:
     https://gitlab.com/libeigen/eigen/-/issues/2000
  */

  em_t a;
  a.setConstant(1);
  const auto ptr = a.data();

  w_t b (std::move(a));
  for (auto j=0; j<4; ++j){
    ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 6.);
  }
  ASSERT_TRUE( b.data()->data() != ptr );
  // ASSERT_TRUE( a.data() == nullptr );
  ASSERT_TRUE( b.data() != nullptr );
}


TEST(containers_dense_matrix_eigen_static, CopyConstructor)
{
  w_t a;
  a.data()->setConstant(1.);

  w_t b(a);
  for (auto j=0; j<a.extent(1); ++j){
    ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 6.);
  }

  // change b should not affect a
  b.data()->setConstant(2.);
  for (auto j=0; j<a.extent(1); ++j){
    ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 12.);
    ASSERT_EQ(a.data()->col(j).lpNorm<1>(), 6.);
  }

  ASSERT_TRUE( a.data()->data() != b.data()->data() );
  ASSERT_TRUE( a.data()->data() != nullptr );
  ASSERT_TRUE( b.data()->data() != nullptr );
}

TEST(containers_dense_matrix_eigen_static, MoveConstructor)
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
  for (auto j=0; j<4; ++j){
    ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 6.);
  }

  // ASSERT_TRUE( a.data()->data() == nullptr );
  ASSERT_TRUE( b.data()->data() != ptr );
}

TEST(containers_dense_matrix_eigen_static, MoveAssign)
{
  /* move semantics for static do not invalidate
     the moved-from object, they just copy the data, see:
     https://gitlab.com/libeigen/eigen/-/issues/2000
  */

  w_t b;
  b.data()->setConstant(1.);
  auto tmp = b.data()->data();
  for (auto j=0; j<4; ++j){
    ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 6.);
  }
  {
    w_t a;
    tmp = a.data()->data();
    a.data()->setConstant(2.);
    ASSERT_TRUE( tmp != nullptr );
    b = std::move(a);
    for (auto j=0; j<4; ++j){
      ASSERT_EQ(b.data()->col(j).lpNorm<1>(), 12.);
    }

    ASSERT_TRUE( b.data()->data() != tmp );

    // waiting for: https://gitlab.com/libeigen/eigen/-/issues/2000
    //ASSERT_TRUE( a.data()->data() == nullptr );
  }
  ASSERT_TRUE( b.data()->data() != tmp );
  ASSERT_TRUE( b.data()->data() != nullptr );
}
