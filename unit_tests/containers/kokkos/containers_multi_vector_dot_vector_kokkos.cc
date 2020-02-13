
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"

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
  static_assert( ::pressio::containers::details::traits<wv>::is_static == 0, "" );
  static_assert( !::pressio::containers::meta::is_multi_vector_wrapper_kokkos<wv>::value, "" );
  static_assert( ::pressio::containers::meta::is_vector_wrapper_kokkos<wv>::value, "" );

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
    // Construct Vector
    // ------------------
    // create host
    k1d_h b_h("bh", Nr);
    // fill host matrix
    b_h(0) = 1.;
    b_h(1) = 1.;
    b_h(2) = 1.;
    b_h(3) = 1.;
    b_h(4) = 2.;
    b_h(5) = 1.;

    // create device
    k1d_d b_d("bd", Nr);
    // copy host -> device
    Kokkos::deep_copy(b_d, b_h);
    // wrap device
    wv vb_d(b_d);

    // ---- do operation ---

    auto c_d = ::pressio::containers::ops::dot(mvA_d, vb_d);
    using expected_ret_t = ::pressio::containers::Vector<k1d_d>;
    static_assert( std::is_same< decltype(c_d), expected_ret_t>::value, "");

    // create host mirror of C_d
    using c_h_t = typename ::pressio::containers::details::traits<wv>::host_mirror_t;
    c_h_t c_h("ch", Nc);
    Kokkos::deep_copy(c_h, *c_d.data());

    EXPECT_EQ(c_h.extent(0), Nc);
    EXPECT_EQ(c_d.data()->extent(0), Nc);
    EXPECT_DOUBLE_EQ( c_h(0), 6.);
    EXPECT_DOUBLE_EQ( c_h(1), 6.);
    EXPECT_DOUBLE_EQ( c_h(2), 6.);
  }//end method
};

TEST(containers_multi_vector_kokkos, dot_vector_dLL){
  // kokkos initialize and finalize already set from environment, see CMakeList

  using d_lay = Kokkos::LayoutLeft;
  RunTest<d_lay> r1;
}

TEST(containers_multi_vector_kokkos, dot_vector_dLR){
  // kokkos initialize and finalize already set from environment, see CMakeList

  using d_lay = Kokkos::LayoutRight;
  RunTest<d_lay> r1;
}
