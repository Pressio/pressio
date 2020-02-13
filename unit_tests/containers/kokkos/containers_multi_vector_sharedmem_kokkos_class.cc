
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"

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
    a[i] = 2.0*i;
  }
};

TEST(containers_multi_vector_sharedmem_kokkos_class, Constructor)
{
  using namespace pressio;

  // kokkos initialize and finalize already set from environment, see CMakeList

  // number of elements
  const int N = 10;

  using view_type = Kokkos::View<double**>;
  view_type a ("A", N, 3);
  //Kokkos::parallel_for ("HelloWorld",15, hello_world());

  using mvec_t = containers::MultiVector<view_type>;
  static_assert( containers::details::traits<mvec_t>::is_static == 0, "" );
  static_assert( containers::meta::is_multi_vector_wrapper_kokkos<mvec_t>::value, "" );
  static_assert( !containers::meta::is_vector_wrapper_kokkos<mvec_t>::value, "" );

  mvec_t g(a);
}


TEST(containers_multi_vector_sharedmem_kokkos, staticCheck)
{
  using namespace pressio;

  {
  using v2_t = Kokkos::View<double**>;
  using mvec_t = containers::MultiVector<v2_t>;
  static_assert( containers::meta::is_multi_vector_wrapper_kokkos<mvec_t>::value, "" );
  static_assert( !containers::meta::is_vector_wrapper_kokkos<mvec_t>::value, "" );
  }

  using v3_t = Kokkos::View<double***>;
  using mvec_t = containers::MultiVector<v3_t>;
  static_assert( !containers::meta::is_multi_vector_wrapper_kokkos<mvec_t>::value, "" );
  static_assert( !containers::meta::is_vector_wrapper_kokkos<mvec_t>::value, "" );
}
