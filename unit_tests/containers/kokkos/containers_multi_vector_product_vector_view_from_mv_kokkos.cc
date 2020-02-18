
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

template <typename d_layout>
struct RunTest{
  // number of elements
  const int Nr = 6;
  const int Nc = 3;

  using exe_space = Kokkos::DefaultExecutionSpace;

  // MV typedefs
  using k2d_d = Kokkos::View<double**, d_layout, exe_space>;
  using k2d_h = typename k2d_d::HostMirror;

  using wmv = ::pressio::containers::MultiVector<k2d_d>;
  static_assert( ::pressio::containers::details::traits<wmv>::is_static == 0, "" );
  static_assert( ::pressio::containers::meta::is_multi_vector_wrapper_kokkos<wmv>::value, "" );
  static_assert( !::pressio::containers::meta::is_vector_wrapper_kokkos<wmv>::value, "" );

  // vector typedefs
  using k1d_d = Kokkos::View<double*, d_layout, exe_space>;
  using k1d_h = typename k1d_d::HostMirror;
  using wv = ::pressio::containers::Vector<k1d_d>;

  RunTest(){
    // ------------------
    // Construct MV data
    // ------------------
    // create host matrix
    k2d_h A_h("Ah", Nr, Nc);
    // fill host matrix
    A_h(0,0) = 1.; A_h(0,1) = 2.; A_h(0,2) = 3.;
    A_h(1,0) = 3.; A_h(1,1) = 2.; A_h(1,2) = 1.;
    A_h(2,0) = 0.; A_h(2,1) = 0.; A_h(2,2) = 1.;
    A_h(3,0) = 0.; A_h(3,1) = 1.; A_h(3,2) = 0.;
    A_h(4,0) = 1.; A_h(4,1) = 0.; A_h(4,2) = 0.;
    A_h(5,0) = 0.; A_h(5,1) = 1.; A_h(5,2) = 1.;

    // create device matrix
    k2d_d A_d("Ad", Nr, Nc);
    // copy host -> device
    Kokkos::deep_copy(A_d, A_h);

    // wrap device matrix
    wmv mvA_d(A_d);

    // ------------------
    // Construct MV2
    // ------------------
    // create host
    k2d_h B_h("bh", Nc, 2);
    // fill host
    B_h(0,0) = 1.; B_h(0,1) = 2.;
    B_h(1,0) = 1.; B_h(1,1) = 2.;
    B_h(2,0) = 1.; B_h(2,1) = 2.;

    // create device
    k2d_d B_d("bd", Nc, 2);
    // copy host -> device
    Kokkos::deep_copy(B_d, B_h);
    // wrap device
    wmv mvB_d(B_d);

    // ---- do operation ---

    {
      // get a view of the 0-th col
      const auto colView = pressio::containers::viewColumnVector(mvB_d, 0);

      auto c_d = ::pressio::containers::ops::product(mvA_d, colView);
      using expected_ret_t = ::pressio::containers::Vector<k1d_d>;
      static_assert( std::is_same< decltype(c_d), expected_ret_t>::value, "");

      // create host mirror of C_d
      using c_h_t = typename ::pressio::containers::details::traits<wv>::host_mirror_t;
      c_h_t c_h("ch", Nr);
      Kokkos::deep_copy(c_h, *c_d.data());

      EXPECT_EQ(c_h.extent(0), Nr);
      EXPECT_EQ(c_d.data()->extent(0), Nr);
      EXPECT_DOUBLE_EQ( c_h(0), 6.);
      EXPECT_DOUBLE_EQ( c_h(1), 6.);
      EXPECT_DOUBLE_EQ( c_h(2), 1.);
      EXPECT_DOUBLE_EQ( c_h(3), 1.);
      EXPECT_DOUBLE_EQ( c_h(4), 1.);
      EXPECT_DOUBLE_EQ( c_h(5), 2.);
    }


    {
      // get a view of the 1-st col
      const auto colView = pressio::containers::viewColumnVector(mvB_d, 1);

      auto c_d = ::pressio::containers::ops::product(mvA_d, colView);
      using expected_ret_t = ::pressio::containers::Vector<k1d_d>;
      static_assert( std::is_same< decltype(c_d), expected_ret_t>::value, "");

      // create host mirror of C_d
      using c_h_t = typename ::pressio::containers::details::traits<wv>::host_mirror_t;
      c_h_t c_h("ch", Nr);
      Kokkos::deep_copy(c_h, *c_d.data());

      EXPECT_EQ(c_h.extent(0), Nr);
      EXPECT_EQ(c_d.data()->extent(0), Nr);
      EXPECT_DOUBLE_EQ( c_h(0), 12.);
      EXPECT_DOUBLE_EQ( c_h(1), 12.);
      EXPECT_DOUBLE_EQ( c_h(2), 2.);
      EXPECT_DOUBLE_EQ( c_h(3), 2.);
      EXPECT_DOUBLE_EQ( c_h(4), 2.);
      EXPECT_DOUBLE_EQ( c_h(5), 4.);
    }

  }//end method
};

TEST(containers_multi_vector_kokkos, dot_vector_view_from_mv_dLL){
  // kokkos initialize and finalize already set from environment

  using d_lay = Kokkos::LayoutLeft;
  RunTest<d_lay> r1;
}
