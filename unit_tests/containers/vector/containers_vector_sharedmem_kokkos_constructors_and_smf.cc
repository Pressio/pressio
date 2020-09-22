
#include <gtest/gtest.h>
#include "pressio_containers.hpp"
#include <KokkosBlas1_fill.hpp>
#include "KokkosBlas1_nrm1.hpp"

using view_t = Kokkos::View<double*, Kokkos::HostSpace>;
using w_t = pressio::containers::Vector<view_t>;

TEST(containers_vector_sharedmem_kokkos, Constructor1)
{
  view_t a ("A", 15);
  KokkosBlas::fill( a, 1. );
  w_t b(a);
  ASSERT_EQ( KokkosBlas::nrm1(*b.data()), 15.);
  // pressio uses copy semantics so the wrapped should
  // NOT be the same as the move-from
  ASSERT_TRUE( a.data() != b.data()->data() );
}

TEST(containers_vector_sharedmem_kokkos, Constructor2)
{
  view_t a ("A", 15);
  KokkosBlas::fill( a, 1. );
  ASSERT_EQ( a.use_count(), 1 );
  w_t b(std::move(a));
  ASSERT_EQ( KokkosBlas::nrm1(*b.data()), 15.);
  ASSERT_EQ( a.use_count(), 0 );
  ASSERT_TRUE( a.data() == b.data()->data() );
 }

TEST(containers_vector_sharedmem_kokkos, Constructor3)
{
  w_t b("label", 15);
  ASSERT_EQ( b.extent(0), 15 );
}

TEST(containers_vector_sharedmem_kokkos, Constructor4)
{
  w_t b(15);
  ASSERT_EQ( b.extent(0), 15 );
}

TEST(containers_vector_sharedmem_kokkos, CopyConstructor)
{
  w_t a(15);
  KokkosBlas::fill( *a.data(), 1. );
  w_t b(a);
  ASSERT_EQ( KokkosBlas::nrm1(*b.data()), 15.);
  ASSERT_TRUE( a.data()->data() != b.data()->data() );
}

TEST(containers_vector_sharedmem_kokkos, CopyAssign)
{
  w_t a(15);
  KokkosBlas::fill( *a.data(), 1. );
  w_t b(15);
  b = a;
  ASSERT_EQ( KokkosBlas::nrm1(*b.data()), 15.);
  ASSERT_TRUE( a.data()->data() != b.data()->data() );
}

// TEST(containers_vector_sharedmem_kokkos, MoveConstructor)
// {
//   w_t a("a",3);
//   ASSERT_EQ( a.data()->use_count(), 1 );
//   w_t b(std::move(a));
//   // ASSERT_TRUE( a.data() == nullptr );
//   ASSERT_EQ( a.data()->use_count(), 0 );
//   ASSERT_TRUE( a.data()->data() == b.data()->data() );
// }

// TEST(containers_vector_sharedmem_kokkos, MoveAssign)
// {
//   view_t a("a",3);
//   ASSERT_EQ( a.use_count(), 1 );
//   view_t b("b",3);
//   b = std::move(a);
//   ASSERT_EQ( a.use_count(), 0 );
//   //ASSERT_TRUE( a.data()->data() != b.data()->data() );
// }
