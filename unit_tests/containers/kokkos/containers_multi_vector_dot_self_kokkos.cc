
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

template <typename d_layout>
struct RunTest{
  // number of elements
  const int Nr = 6;
  const int Nc = 4;

  using exe_space = Kokkos::DefaultExecutionSpace;
  using k2d_d = Kokkos::View<double**, d_layout, exe_space>;
  using k2d_h = typename k2d_d::HostMirror;

  using wmv = ::pressio::containers::MultiVector<k2d_d>;
  static_assert( ::pressio::containers::details::traits<wmv>::is_static == 0, "" );
  static_assert( ::pressio::containers::meta::is_multi_vector_wrapper_kokkos<wmv>::value, "" );
  static_assert( !::pressio::containers::meta::is_vector_wrapper_kokkos<wmv>::value, "" );

  RunTest(){
    // create host matrix
    k2d_h A_h("Ah", Nr, Nc);
    // fill host matrix
    A_h(0,0)=0.; A_h(0,1)=1.; A_h(0,2)=2.; A_h(0,3)=3.;
    A_h(1,0)=1.; A_h(1,1)=0.; A_h(1,2)=1.; A_h(1,3)=2.;
    A_h(2,0)=0.; A_h(2,1)=2.; A_h(2,2)=0.; A_h(2,3)=1.;
    A_h(3,0)=0.; A_h(3,1)=2.; A_h(3,2)=2.; A_h(3,3)=1.;
    A_h(4,0)=0.; A_h(4,1)=3.; A_h(4,2)=0.; A_h(4,3)=0.;
    A_h(5,0)=4.; A_h(5,1)=0.; A_h(5,2)=1.; A_h(5,3)=0.;

    // create device matrix
    k2d_d A_d("Ad", Nr, Nc);
    // copy host -> device
    Kokkos::deep_copy(A_d, A_h);

    // wrap device matrix
    wmv mvA_d(A_d);

    using ret_t = ::pressio::containers::Matrix<k2d_d>;
    // self-dot
    constexpr auto beta  = ::pressio::utils::constants::zero<sc_t>();
    constexpr auto alpha = ::pressio::utils::constants::one<sc_t>();
    auto C_d = ::pressio::containers::ops::product<ret_t>(::pressio::transpose(), ::pressio::nontranspose(), alpha, MV);

    using expected_ret_t = ::pressio::containers::Matrix<k2d_d>;
    static_assert( std::is_same< decltype(C_d), expected_ret_t>::value, "");

    // create host mirror of C_d
    using C_h_t = typename ::pressio::containers::details::traits<wmv>::host_mirror_t;
    C_h_t C_h("Ch", Nc, Nc);
    Kokkos::deep_copy(C_h, *C_d.data());

    // -- check solution --
    Eigen::MatrixXd TT(4,4);
    TT(0,0)=17.; TT(0,1)=0.; TT(0,2)=5.; TT(0,3)=2.;
    TT(1,0)=0.; TT(1,1)=18.; TT(1,2)=6.; TT(1,3)=7.;
    TT(2,0)=5.; TT(2,1)=6.; TT(2,2)=10.; TT(2,3)=10.;
    TT(3,0)=2.; TT(3,1)=7.; TT(3,2)=10.; TT(3,3)=15.;
    for (auto i=0; i<Nc; i++){
      for (auto j=0; j<Nc; j++){
	EXPECT_NEAR( TT(i,j), C_h(i,j), 1e-12);
      }
    }
  }
};

TEST(containers_multi_vector_kokkos, dotSelf_dLL){
  // kokkos initialize and finalize already set from environment, see CMakeList

  using d_lay = Kokkos::LayoutLeft;
  RunTest<d_lay> r1;
}

TEST(containers_multi_vector_kokkos, dotSelf_dLR){
  // kokkos initialize and finalize already set from environment, see CMakeList

  using d_lay = Kokkos::LayoutRight;
  RunTest<d_lay> r1;
}
