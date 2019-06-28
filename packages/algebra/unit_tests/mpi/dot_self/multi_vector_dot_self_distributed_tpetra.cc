
#include <gtest/gtest.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_CommHelpers.hpp>
#include "ALGEBRA_ALL"
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

    result is symmetric
    A^T A =  17  0    5   2
	      0  18   6   7
	      5  6   10  10
	      2	 7   10  15
   */

TEST(algebra_ops, TpetraMultiVectorDotSelf){
  using namespace rompp;

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
  using mvec_t = algebra::MultiVector<nat_mvec_t>;
  using mv_device_t = typename algebra::details::traits<mvec_t>::device_t;

  mvec_t A( *mvA ); A.setZero();

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

  using eig_mat = Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic>;
  eig_mat TT(4,4);
  TT(0,0)=17.; TT(0,1)=0.; TT(0,2)=5.; TT(0,3)=2.;
  TT(1,0)=0.; TT(1,1)=18.; TT(1,2)=6.; TT(1,3)=7.;
  TT(2,0)=5.; TT(2,1)=6.; TT(2,2)=10.; TT(2,3)=10.;
  TT(3,0)=2.; TT(3,1)=7.; TT(3,2)=10.; TT(3,3)=15.;

  auto C = algebra::ops::dot_self(A);
  for (auto i=0; i<4; i++){
    for (auto j=0; j<2; j++){
      EXPECT_NEAR( TT(i,j), C(i,j), 1e-12);
    }
  }

  algebra::Matrix<eig_mat> C2(4,4);
  algebra::ops::dot_self(A, C2);
  for (auto i=0; i<4; i++){
    for (auto j=0; j<4; j++){
      EXPECT_NEAR( TT(i,j), C2(i,j), 1e-12);
    }
  }
}
