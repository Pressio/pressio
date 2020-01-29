
#include "tpetra_only_fixtures.hpp"

TEST_F(tpetraMultiVectorGlobSize9Fixture,
       MVecTpetraProductKokkosVector){
  using namespace pressio;

  using nat_mv_t = typename tpetraMultiVectorGlobSize9Fixture::mvec_t;
  using mvec_t = containers::MultiVector<nat_mv_t>;
  STATIC_ASSERT_IS_CONTAINERS_MULTI_VECTOR_WRAPPER(mvec_t);
  using device_t = typename containers::details::traits<mvec_t>::device_t;

  // construct multivector wrapper
  mvec_t MV( *x_ );
  ::pressio::containers::ops::set_zero(MV);

  // get trilinos tpetra multivector object
  auto trilD = MV.data();
  trilD->sync<Kokkos::HostSpace>();

  /*--------------------------------------------
   * (1): modify the host view and then sync
  //--------------------------------------------*/
  auto v2d = trilD->getLocalView<Kokkos::HostSpace>();
  auto c0 = Kokkos::subview(v2d, Kokkos::ALL(), 0);
  //we are going to change the host view
  trilD->modify<Kokkos::HostSpace>();

  if(rank_==0){
    v2d(0,0) = 3.2; v2d(0,1) = 1.2;
    v2d(1,0) = 1.2;
    v2d(2,1) = 4; v2d(2,3) = 1;
  }
  if(rank_==1){
    v2d(2,2) = 3;
  }
  if(rank_==2){
    v2d(1,1) = 2.2;
  }
  // sync from host to device
  trilD->sync<device_t> ();

  // MyPID    GID
  //     0     0      3.2     1.2     0      0
  //     0     1      1.2      0      0      0
  //     0     2       0       4      0      1
  //     1     3       0       0      0      0
  //     1     4       0       0      0      0
  //     1     5       0       0      3      0
  //     2     6       0       0      0      0
  //     2     7       0       2.2    0      0
  //     2     8       0       0      0      0

  EXPECT_EQ( MV.globalNumVectors(), 4 );
  EXPECT_EQ( MV.localNumVectors(), 4 );
  EXPECT_EQ( MV.globalLength(), 9 );
  EXPECT_EQ( MV.localLength(), 3);
  //----------

  // create the vector to do the product
  // the vector is a wrapper of a kokkos vector
  {
    using k1d_d = Kokkos::View<double*, device_t>;
    using k1d_h = typename k1d_d::HostMirror;
    using vec_t = pressio::containers::Vector<k1d_d>;
    STATIC_ASSERT_IS_CONTAINERS_VECTOR_WRAPPER(vec_t);

    // create host vector
    k1d_h b_h("bh", 4);
    b_h(0) = 1.; b_h(1) = 2.; b_h(2) = 3.; b_h(3) = 4.;
    // create device
    k1d_d b_d("bd", 4);
    // copy host -> device
    Kokkos::deep_copy(b_d, b_h);
    // wrap device vector
    vec_t b(b_d);
    // do product, returning a tpetra vector
    auto res = containers::ops::product(MV, b);
    auto cc2 = res.data()->getLocalView<Kokkos::HostSpace>();
    auto cc = Kokkos::subview (cc2, Kokkos::ALL (), 0);
    if (rank_==0){
      EXPECT_DOUBLE_EQ( cc(0), 5.6);
      EXPECT_DOUBLE_EQ( cc(1), 1.2);
      EXPECT_DOUBLE_EQ( cc(2), 12.0);
    }
    if (rank_==1){
      EXPECT_DOUBLE_EQ( cc(0), 0.0);
      EXPECT_DOUBLE_EQ( cc(1), 0.0);
      EXPECT_DOUBLE_EQ( cc(2), 9.0);
    }
    if (rank_==2){
      EXPECT_DOUBLE_EQ( cc(0), 0.);
      EXPECT_DOUBLE_EQ( cc(1), 4.4);
      EXPECT_DOUBLE_EQ( cc(2), 0.);
    }
  }
}
