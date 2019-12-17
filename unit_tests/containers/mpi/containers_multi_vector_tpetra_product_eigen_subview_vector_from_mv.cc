
#include "tpetra_only_fixtures.hpp"

TEST_F(tpetraMultiVectorGlobSize9Fixture,
       MVecTpetraProductEigenVectorViewFromMV){
  using namespace pressio;

  // aliases
  using nat_mv_t = typename tpetraMultiVectorGlobSize9Fixture::mvec_t;
  using mvec_t = containers::MultiVector<nat_mv_t>;
  STATIC_ASSERT_IS_CONTAINERS_MULTI_VECTOR_WRAPPER(mvec_t);
  using device_t = typename containers::details::traits<mvec_t>::device_t;

  // construct multivector wrapper
  mvec_t MV( *x_ );
  MV.setZero();
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
    v2d(2,1) = 4;
  }
  if(rank_==1){
    v2d(2,2) = 3;
  }
  // sync from host to device
  trilD->sync<device_t> ();

// MyPID    GID
//     0     0      3.2     1.2     0      0
//     0     1      1.2      0      0      0
//     0     2       0       4      0      0
//     1     3       0       0      0      0
//     1     4       0       0      0      0
//     1     5       0       0      3      0
//     2     6       0       0      0      0
//     2     7       0       0      0      0
//     2     8       0       0      0      0

  EXPECT_EQ( MV.globalNumVectors(), 4 );
  EXPECT_EQ( MV.localNumVectors(), 4 );
  EXPECT_EQ( MV.globalLength(), 9 );
  EXPECT_EQ( MV.localLength(), 3);
  //----------

  using eig_mat_t = Eigen::MatrixXd;
  eig_mat_t bn(4,2);
  using mv_t = containers::MultiVector<eig_mat_t>;
  mv_t B(bn);
  B(0, 0) = 1.; B(0, 1) = 2.;
  B(1, 0) = 2.; B(1, 1) = 4.;
  B(2, 0) = 3.; B(2, 1) = 6.;
  B(3, 0) = 4.; B(3, 1) = 8.;

  {
    const auto colVec = pressio::containers::viewColumnVector(B, 0);
    auto res = containers::ops::product(MV, colVec);
    auto cc2 = res.data()->getLocalView<Kokkos::HostSpace>();
    auto cc = Kokkos::subview (cc2, Kokkos::ALL (), 0);
    if (rank_==0){
      EXPECT_DOUBLE_EQ( cc(0), 5.6);
      EXPECT_DOUBLE_EQ( cc(1), 1.2);
      EXPECT_DOUBLE_EQ( cc(2), 8.0);
    }
    if (rank_==1){
      EXPECT_DOUBLE_EQ( cc(0), 0.0);
      EXPECT_DOUBLE_EQ( cc(1), 0.0);
      EXPECT_DOUBLE_EQ( cc(2), 9.0);
    }
    if (rank_==2){
      EXPECT_DOUBLE_EQ( cc(0), 0.);
      EXPECT_DOUBLE_EQ( cc(1), 0.);
      EXPECT_DOUBLE_EQ( cc(2), 0.);
    }

    res.data()->putScalar(0.);
    containers::ops::product(MV, colVec, res);
    {
      auto cc2 = res.data()->getLocalView<Kokkos::HostSpace>();
      auto cc = Kokkos::subview (cc2, Kokkos::ALL (), 0);
      if (rank_==0){
	EXPECT_DOUBLE_EQ( cc(0), 5.6);
	EXPECT_DOUBLE_EQ( cc(1), 1.2);
	EXPECT_DOUBLE_EQ( cc(2), 8.0);
      }
      if (rank_==1){
	EXPECT_DOUBLE_EQ( cc(0), 0.0);
	EXPECT_DOUBLE_EQ( cc(1), 0.0);
	EXPECT_DOUBLE_EQ( cc(2), 9.0);
      }
      if (rank_==2){
	EXPECT_DOUBLE_EQ( cc(0), 0.);
	EXPECT_DOUBLE_EQ( cc(1), 0.);
	EXPECT_DOUBLE_EQ( cc(2), 0.);
      }
    }
  }


  {
    const auto colVec = pressio::containers::viewColumnVector(B, 1);
    auto res = containers::ops::product(MV, colVec);
    auto cc2 = res.data()->getLocalView<Kokkos::HostSpace>();
    auto cc = Kokkos::subview (cc2, Kokkos::ALL (), 0);
    if (rank_==0){
      EXPECT_DOUBLE_EQ( cc(0), 11.2);
      EXPECT_DOUBLE_EQ( cc(1), 2.4);
      EXPECT_DOUBLE_EQ( cc(2), 16.0);
    }
    if (rank_==1){
      EXPECT_DOUBLE_EQ( cc(0), 0.0);
      EXPECT_DOUBLE_EQ( cc(1), 0.0);
      EXPECT_DOUBLE_EQ( cc(2), 18.0);
    }
    if (rank_==2){
      EXPECT_DOUBLE_EQ( cc(0), 0.);
      EXPECT_DOUBLE_EQ( cc(1), 0.);
      EXPECT_DOUBLE_EQ( cc(2), 0.);
    }

    res.data()->putScalar(0.);
    containers::ops::product(MV, colVec, res);
    {
      auto cc2 = res.data()->getLocalView<Kokkos::HostSpace>();
      auto cc = Kokkos::subview (cc2, Kokkos::ALL (), 0);
      if (rank_==0){
	EXPECT_DOUBLE_EQ( cc(0), 11.2);
	EXPECT_DOUBLE_EQ( cc(1), 2.4);
	EXPECT_DOUBLE_EQ( cc(2), 16.0);
      }
      if (rank_==1){
	EXPECT_DOUBLE_EQ( cc(0), 0.0);
	EXPECT_DOUBLE_EQ( cc(1), 0.0);
	EXPECT_DOUBLE_EQ( cc(2), 18.0);
      }
      if (rank_==2){
	EXPECT_DOUBLE_EQ( cc(0), 0.);
	EXPECT_DOUBLE_EQ( cc(1), 0.);
	EXPECT_DOUBLE_EQ( cc(2), 0.);
      }
    }
  }
}
