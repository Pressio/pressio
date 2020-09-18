
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

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
