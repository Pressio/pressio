
#include "tpetra_only_fixtures.hpp"

TEST_F(tpetraMultiVectorGlobSize9Fixture,
       MVecTpetraProductEigenVector){
  using namespace pressio;

  // aliases
  using nat_mv_t = typename tpetraMultiVectorGlobSize9Fixture::mvec_t;
  using mvec_t = containers::MultiVector<nat_mv_t>;
  STATIC_ASSERT_IS_CONTAINERS_MULTI_VECTOR_WRAPPER(mvec_t);
  using device_t = typename containers::details::traits<mvec_t>::device_t;
  // using device_mem_space = typename device_t::memory_space;
  // using mag_type = typename containers::details::traits<mvec_t>::mag_t;

  // construct multivector wrapper
  mvec_t MV( *x_ );
  ::pressio::ops::set_zero(MV);
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

  EXPECT_EQ( MV.numVectors(), 4 );
  EXPECT_EQ( MV.numVectorsLocal(), 4 );
  EXPECT_EQ( MV.extent(0), 9 );
  EXPECT_EQ( MV.extentLocal(0), 3);
  //----------

  using eigv_t = Eigen::Matrix<double,4,1>;
  eigv_t bn;
  using vec_t = containers::Vector<eigv_t>;
  STATIC_ASSERT_IS_CONTAINERS_VECTOR_WRAPPER(vec_t);
  vec_t b(bn);
  b[0] = 1.;
  b[1] = 2.;
  b[2] = 3.;
  b[3] = 4.;

  auto res = ops::product(MV, b);
  auto cc2 = res.data()->getLocalView<Kokkos::HostSpace>();
  auto cc = Kokkos::subview (cc2, Kokkos::ALL (), 0);
  // if(rank_==0){
  //   for (auto i=0;i<res.localSize(); i++)
  //     std::cout << cc(i) << std::endl;
  // }
  //res.data()->Print(std::cout);

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
