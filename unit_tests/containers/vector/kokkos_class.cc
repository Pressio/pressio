
#include <gtest/gtest.h>
#include "pressio_containers.hpp"
#include <KokkosBlas1_fill.hpp>
#include "KokkosBlas1_nrm1.hpp"

template <typename T>
struct InitView {
  T a;

  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  InitView (T a_) :
    a (a_){}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const int i) const {
    // Acesss the View just like a Fortran array.  The layout depends
    // on the View's memory space, so don't rely on the View's
    // physical memory layout unless you know what you're doing.
    a(i) = 2.0*i;
  }
};

TEST(containers_vector_sharedmem_kokkos, dynamicWrapper)
{
  using namespace pressio;
  using view_type = Kokkos::View<double*>;
  using myvec_t = containers::Vector<view_type>;
  static_assert(containers::predicates::is_dynamic_vector_wrapper_kokkos<myvec_t>::value, "" );
  static_assert(!containers::predicates::is_static_vector_wrapper_kokkos<myvec_t>::value, "" );
  static_assert(containers::details::traits<myvec_t>::is_static == 0, "" );
}

TEST(containers_vector_sharedmem_kokkos, staticWrapper)
{
  using namespace pressio;
  using view_type = Kokkos::View<double[3]>;
  using myvec_t = containers::Vector<view_type>;
  static_assert(!containers::predicates::is_dynamic_vector_wrapper_kokkos<myvec_t>::value, "" );
  static_assert(containers::predicates::is_static_vector_wrapper_kokkos<myvec_t>::value, "" );
  static_assert(containers::details::traits<myvec_t>::is_static == 1, "" );
}

TEST(containers_vector_sharedmem_kokkos_class, Constructor)
{
  using namespace pressio;

  // kokkos initialize and finalize already set from environment, see CMakeList
  const int N = 10;

  using view_type = Kokkos::View<double*>;
  view_type a ("A", N);
  //Kokkos::parallel_for ("HelloWorld",15, hello_world());

  using myvec_t = containers::Vector<view_type>;
  static_assert( containers::details::traits<myvec_t>::is_static == 0, "" );
  static_assert( containers::predicates::is_vector_wrapper_kokkos<myvec_t>::value, "" );
  static_assert( !containers::predicates::is_multi_vector_wrapper_kokkos<myvec_t>::value, "" );

  myvec_t g(a);

  Kokkos::parallel_for (N, InitView<view_type>( *g.data() ));
  // double sum = 0;
  // Kokkos::parallel_reduce (N, ReduceFunctor<view_type>(a), sum);
  // printf ("Result: %f\n", sum);

  using view_type2 = Kokkos::View<double[3]>;
  using myvec_t2 = containers::Vector<view_type2>;
  static_assert( containers::details::traits<myvec_t2>::is_static == 1, "" );

}



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

  ASSERT_EQ( b.data()->use_count(), 1 );
  ASSERT_TRUE( a.data() == b.data()->data() );

  // just for testing, dangerous to do in reality
  a.~view_t();
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
  // when we copyConstruct a kokkos wrapper
  // we want to make sure that the new object does not
  // view the previous one
  w_t a(15);
  KokkosBlas::fill( *a.data(), 1. );
  ASSERT_TRUE( a.data()->use_count() == 1);

  w_t b(a);
  ASSERT_EQ( KokkosBlas::nrm1(*b.data()), 15.);
  ASSERT_TRUE( a.data()->data() != b.data()->data() );
  ASSERT_TRUE( a.data()->use_count() == 1);
  ASSERT_TRUE( b.data()->use_count() == 1);

  // change b, which  should NOT change a
  KokkosBlas::fill( *b.data(), 3. );
  ASSERT_EQ( KokkosBlas::nrm1(*b.data()), 45.);
  ASSERT_EQ( KokkosBlas::nrm1(*a.data()), 15.);
}

TEST(containers_vector_sharedmem_kokkos, MoveConstructor)
{
  w_t a("a",15);
  KokkosBlas::fill( *a.data(), 1. );
  ASSERT_EQ( KokkosBlas::nrm1(*a.data()), 15.);
  auto aptr = a.data()->data();

  w_t b(std::move(a));
  auto bptr = b.data()->data();
  ASSERT_EQ( KokkosBlas::nrm1(*b.data()), 15.);
  // if we move, the ptr should be same of a
  ASSERT_TRUE( bptr == aptr );

  // explicitly call destructor JUST for testing purposes
  // in this case this should work because of the move
  a.~w_t();

  ASSERT_TRUE( b.data()->use_count() == 1 );
  ASSERT_EQ( KokkosBlas::nrm1(*b.data()), 15.);
}


TEST(containers_vector_sharedmem_kokkos, MoveAssign)
{
  w_t b("a",15);
  KokkosBlas::fill( *b.data(), 1. );
  ASSERT_EQ( KokkosBlas::nrm1(*b.data()), 15.);
  auto bptr1 = b.data()->data();
  {
    w_t a("a",15);
    KokkosBlas::fill( *a.data(), 2. );
    auto aptr = a.data()->data();

    b = std::move(a);
    auto bptr2 = b.data()->data();
    // if we move, the ptr should be same of a
    ASSERT_TRUE(  bptr2 == aptr );
    ASSERT_FALSE( bptr1 == aptr );
  }

  ASSERT_TRUE( b.data()->use_count() == 1 );
  ASSERT_EQ( KokkosBlas::nrm1(*b.data()), 30.);
}


  // view_t a ("A", 15);
  // std::cout << a.use_count() << std::endl;
  // view_t b (a);
  // std::cout << a.use_count() << std::endl;
  // std::cout << b.use_count() << std::endl;

  // view_t c (std::move(a));
  // std::cout << a.use_count() << std::endl;
  // std::cout << c.use_count() << std::endl;

  // view_t d ("d",15);
  // d = std::move(c);
  // std::cout << a.use_count() << std::endl;
  // std::cout << c.use_count() << std::endl;
  // std::cout << d.use_count() << std::endl;
