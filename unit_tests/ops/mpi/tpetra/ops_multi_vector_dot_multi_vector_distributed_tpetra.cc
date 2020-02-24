
#include <gtest/gtest.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_CommHelpers.hpp>
#include "pressio_containers.hpp"
#include <Tpetra_Map_decl.hpp>

  /*
    MultiVector A
    0 1 2 3
    1 0 1 2
    -------
    0 2 0 1
    0 2 2 1
    -------
    0 3 0 0
    4 0 1 0

    MultiVector B
    1 2
    0 0
    ----
    1 0
    1 0
    ----
    0 2
    0 1


    A^T B =  0 4
	     5 8
	     4 5
	     5 6
   */


TEST(ops, TpetraMultiVectorDotMultiVector){
  using namespace pressio;

  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using nat_mvec_t = Tpetra::MultiVector<>;

  Teuchos::RCP<const tcomm> comm_ =
	Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  auto rank = comm_->getRank();
  auto numProc_ = comm_->getSize();
  EXPECT_EQ(numProc_,3);

  Teuchos::RCP<const map_t> contigMap_ =
    Teuchos::rcp(new map_t(6, 0, comm_));

  auto mvA = std::make_shared<nat_mvec_t>(contigMap_, 4);
  auto mvB = std::make_shared<nat_mvec_t>(contigMap_, 2);

  using mvec_t = containers::MultiVector<nat_mvec_t>;
  using mv_device_t = typename containers::details::traits<mvec_t>::device_t;

  mvec_t A( *mvA ); ::pressio::ops::set_zero(A);
  mvec_t B( *mvB ); ::pressio::ops::set_zero(B);

  /*------------- fill A -------------- */
  {
    // get trilinos tpetra multivector object
    auto trilD = A.data();
    trilD->sync<Kokkos::HostSpace>();

    /*--------------------------------------------
     * (1): modify the host view and then sync
     //--------------------------------------------*/
    auto v2d = trilD->getLocalView<Kokkos::HostSpace>();
    auto c0 = Kokkos::subview(v2d, Kokkos::ALL(), 0);
    //we are going to change the host view
    trilD->modify<Kokkos::HostSpace>();

    if(rank==0){
      v2d(0,0)=0.; v2d(0,1)=1; v2d(0,2)=2; v2d(0,3)=3;
      v2d(1,0)=1.; v2d(1,1)=0; v2d(1,2)=1; v2d(1,3)=2;
    }
    if(rank==1){
      v2d(0,0)=0.; v2d(0,1)=2; v2d(0,2)=0; v2d(0,3)=1;
      v2d(1,0)=0.; v2d(1,1)=2; v2d(1,2)=2; v2d(1,3)=1;
    }
    if(rank==2){
      v2d(0,0)=0.; v2d(0,1)=3; v2d(0,2)=0; v2d(0,3)=0;
      v2d(1,0)=4.; v2d(1,1)=0; v2d(1,2)=1; v2d(1,3)=0;
    }
    // sync from host to device
    trilD->sync<mv_device_t> ();
  }

  {
    // get trilinos tpetra multivector object
    auto trilD = B.data();
    trilD->sync<Kokkos::HostSpace>();

    /*--------------------------------------------
     * (1): modify the host view and then sync
     //--------------------------------------------*/

    auto v2d = trilD->getLocalView<Kokkos::HostSpace>();
    auto c0 = Kokkos::subview(v2d, Kokkos::ALL(), 0);
    //we are going to change the host view
    trilD->modify<Kokkos::HostSpace>();

    if(rank==0){
      v2d(0,0)=1.; v2d(0,1)=2;
      v2d(1,0)=0.; v2d(1,1)=0;
    }
    if(rank==1){
      v2d(0,0)=1.; v2d(0,1)=0;
      v2d(1,0)=1.; v2d(1,1)=0;
    }
    if(rank==2){
      v2d(0,0)=0.; v2d(0,1)=2;
      v2d(1,0)=0.; v2d(1,1)=1;
    }
    // sync from host to device
    trilD->sync<mv_device_t> ();
  }

  using eig_mat = Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic>;
  using eig_mat_w = containers::Matrix<eig_mat>;

  eig_mat TT(4,2);
  TT(0,0) = 0.0; TT(0,1) = 4.0;
  TT(1,0) = 5.0; TT(1,1) = 8.0;
  TT(2,0) = 4.0; TT(2,1) = 5.0;
  TT(3,0) = 5.0; TT(3,1) = 6.0;


  eig_mat_w C2(4,2);
  constexpr auto beta  = ::pressio::utils::constants::zero<scalar_t>();
  constexpr auto alpha = ::pressio::utils::constants::one<scalar_t>();
  ::pressio::ops::product(::pressio::transpose(), 
    ::pressio::nontranspose(), alpha, A, B, beta, C2);

  for (auto i=0; i<4; i++){
    for (auto j=0; j<2; j++){
      EXPECT_NEAR( TT(i,j), C2(i,j), 1e-12);
    }
  }
}
