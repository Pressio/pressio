
#include <gtest/gtest.h>
#include "CORE_VECTOR"
#include<Kokkos_Core.hpp>
#include<Kokkos_Random.hpp>
#include<KokkosBlas1_axpby.hpp>
//#include<KokkosBlas1_dot.hpp>

TEST(core_vector_sharedmem_kokkos_class, axpy){
  // kokkos initialize and finalize already set from environment, see parent CMakeList
  
  const int N = 10;
  using view_type = Kokkos::View<double*>;
  view_type x ("x", N);
  view_type y ("y", N);
  // create host mirrors
  typename view_type::HostMirror xH = Kokkos::create_mirror(x);
  typename view_type::HostMirror yH = Kokkos::create_mirror(y);

  // // fill device views
  // Kokkos::fill_random(x, Kokkos::impl::rand_pool, 1.1);
  // Kokkos::fill_random(y, Kokkos::impl::rand_pool, 1.1);

  // fill y
  Kokkos::parallel_for (N,
			KOKKOS_LAMBDA (const int i) {
			  y[i] = 1.1;
			  x[i] = 1.0;
			});
  
  // copy to host
  Kokkos::deep_copy (yH, y);
  Kokkos::deep_copy (xH, x);
  
  // y = a*x + y
  double a = 3.0;
  KokkosBlas::axpy(a,x,y);

  Kokkos::deep_copy (yH, y);
  
  for (size_t i=0; i<yH.size(); i++)
    std::cout << yH[i] << ' ';
  
  
  // using myvec_t = core::Vector<view_type>;
  // static_assert( core::details::traits<myvec_t>::is_static == 0, "" );
  // myvec_t g(a);
  
  // using view_type2 = Kokkos::View<double[11]>;
  // using myvec_t2 = core::Vector<view_type2>;
  // static_assert( core::details::traits<myvec_t2>::is_static == 1, "" );  
}
