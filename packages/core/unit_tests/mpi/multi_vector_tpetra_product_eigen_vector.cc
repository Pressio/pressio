
#include "tpetra_only_fixtures.hpp"

TEST_F(tpetraMultiVectorGlobSize9Fixture,
       MVecTpetraProductEigenVector){
  using namespace rompp;

  // aliases
  using nat_mv_t = typename tpetraMultiVectorGlobSize9Fixture::mvec_t;
  using mvec_t = core::MultiVector<nat_mv_t>;
  STATIC_ASSERT_IS_CORE_MULTI_VECTOR_WRAPPER(mvec_t);
  using device_t = typename core::details::traits<mvec_t>::device_t;
  // using device_mem_space = typename device_t::memory_space;
  // using mag_type = typename core::details::traits<mvec_t>::mag_t;

  // construct multivector wrapper
  mvec_t MV( *x_ );
  MV.setZero();
  // get trilinos tpetra multivector object
  auto trilD = MV.data();
  
  /*--------------------------------------------
   * (1): modify the host view and then sync
   * most likely, host and device will be same unless we run CUDA
   * so in theory we should not worry about syncing but it 
   * does not hurt to do it anyway
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

  using eigv_t = Eigen::Matrix<double,4,1>;
  eigv_t bn;
  using vec_t = core::Vector<eigv_t>;
  STATIC_ASSERT_IS_CORE_VECTOR_WRAPPER(vec_t);
  vec_t b(bn);
  b[0] = 1.;
  b[1] = 2.;
  b[2] = 3.;
  b[3] = 4.;

  auto res = core::ops::product(MV, b); 
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



  // // Now compute the inf-norms of the columns of X.
  // Kokkos::DualView<mag_type*, device_t> norms ("norms", 4);
  // norms.template modify<device_t>();
  // trilD->norm1(norms.template view<device_t>());
  // norms.template sync<Kokkos::HostSpace>();
  // // check on host
  // for (size_t k = 0; k < 4; ++k) {
  //   std::cout << norms.h_view(k) << std::endl;
  // }
