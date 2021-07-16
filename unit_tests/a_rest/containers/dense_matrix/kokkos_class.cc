
#include <gtest/gtest.h>
#include "pressio_containers.hpp"
#include <KokkosBlas1_fill.hpp>
#include "KokkosBlas1_nrm1.hpp"

using view1d_t = Kokkos::View<double*, Kokkos::HostSpace>;
using view2d_t = Kokkos::View<double**, Kokkos::HostSpace>;
using w_t = pressio::containers::DenseMatrix<view2d_t>;

TEST(containers_dense_matrix_kokkos, Constructor1)
{
  view2d_t a ("A", 15, 3);
  KokkosBlas::fill( a, 1. );
  w_t b(a);
  view1d_t norms("n",3); KokkosBlas::nrm1(norms, *b.data());
  for (auto i=0; i<3; ++i) ASSERT_EQ( norms(i), 15.);

  // pressio uses copy semantics so the wrapped should
  // NOT be the same as the move-from
  ASSERT_TRUE( a.data() != b.data()->data() );
}

TEST(containers_dense_matrix_kokkos, Constructor2)
{
  view2d_t a ("A", 15, 3);
  KokkosBlas::fill( a, 1. );
  ASSERT_EQ( a.use_count(), 1 );

  w_t b(std::move(a));
  view1d_t norms("n",3); KokkosBlas::nrm1(norms, *b.data());
  for (auto i=0; i<3; ++i) ASSERT_EQ( norms(i), 15.);

  ASSERT_EQ( b.data()->use_count(), 1 );
  ASSERT_TRUE( a.data() == b.data()->data() );

  // just for testing, dangerous to do in reality
  a.~view2d_t();
 }

TEST(containers_dense_matrix_kokkos, Constructor3)
{
  w_t b("label", 15, 3);
  ASSERT_EQ( b.extent(0), 15 );
  ASSERT_EQ( b.extent(1), 3 );
}

TEST(containers_dense_matrix_kokkos, Constructor4)
{
  w_t b(15,3);
  ASSERT_EQ( b.extent(0), 15 );
  ASSERT_EQ( b.extent(1), 3 );
}

TEST(containers_dense_matrix_kokkos, CopyConstructor)
{
  view1d_t norms("n",3); 

  // when we copyConstruct a kokkos wrapper
  // we want to make sure that the new object does not
  // view the previous one
  w_t a(15,3);
  KokkosBlas::fill( *a.data(), 1. );
  ASSERT_TRUE( a.data()->use_count() == 1);

  w_t b(a);
  KokkosBlas::nrm1(norms, *b.data());
  for (auto i=0; i<3; ++i) ASSERT_EQ( norms(i), 15.);

  ASSERT_TRUE( a.data()->data() != b.data()->data() );
  ASSERT_TRUE( a.data()->use_count() == 1);
  ASSERT_TRUE( b.data()->use_count() == 1);

  // change b, which  should NOT change a
  KokkosBlas::fill( *b.data(), 3. );

  KokkosBlas::nrm1(norms, *b.data());
  for (auto i=0; i<3; ++i) ASSERT_EQ( norms(i), 45.);

  KokkosBlas::nrm1(norms, *a.data());
  for (auto i=0; i<3; ++i) ASSERT_EQ( norms(i), 15.);
}

TEST(containers_dense_matrix_kokkos, MoveConstructor)
{
  view1d_t norms("n",3); 

  w_t a("a",15,3);
  KokkosBlas::fill( *a.data(), 1. );
  KokkosBlas::nrm1(norms, *a.data());
  for (auto i=0; i<3; ++i) ASSERT_EQ( norms(i), 15.);
  auto aptr = a.data()->data();

  w_t b(std::move(a));
  auto bptr = b.data()->data();
  KokkosBlas::nrm1(norms, *b.data());
  for (auto i=0; i<3; ++i) ASSERT_EQ( norms(i), 15.);

  // if we move, the ptr should be same of a
  ASSERT_TRUE( bptr == aptr );

  // explicitly call destructor JUST for testing purposes
  // in this case this should work because of the move
  a.~w_t();

  ASSERT_TRUE( b.data()->use_count() == 1 );
  KokkosBlas::nrm1(norms, *b.data());
  for (auto i=0; i<3; ++i) ASSERT_EQ( norms(i), 15.);
}


TEST(containers_dense_matrix_kokkos, MoveAssign)
{
  view1d_t norms("n",3); 

  w_t b("a",15,3);
  KokkosBlas::fill( *b.data(), 1. );
  KokkosBlas::nrm1(norms, *b.data());
  for (auto i=0; i<3; ++i) ASSERT_EQ( norms(i), 15.);

  auto bptr1 = b.data()->data();
  {
    w_t a("a",15,3);
    KokkosBlas::fill( *a.data(), 2. );
    auto aptr = a.data()->data();

    b = std::move(a);
    auto bptr2 = b.data()->data();
    // if we move, the ptr should be same of a
    ASSERT_TRUE(  bptr2 == aptr );
    ASSERT_FALSE( bptr1 == aptr );
  }

  ASSERT_TRUE( b.data()->use_count() == 1 );
  KokkosBlas::nrm1(norms, *b.data());
  for (auto i=0; i<3; ++i) ASSERT_EQ( norms(i), 30.);
}
