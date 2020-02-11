
#ifndef QR_UTEST_FIXTURES_HPP_
#define QR_UTEST_FIXTURES_HPP_

#include <gtest/gtest.h>
#include "CONTAINERS_ALL"
#include "qr_r9c4_gold.hpp"

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "Epetra_MpiComm.h"
#include "Eigen/Dense"
#endif

namespace pressio{ namespace qr{ namespace test{
  static constexpr int numVectors_ = 4;
  static constexpr int numRows_ = 9;
}}}// end namespace pressio::qr::test


struct eigenDenseR9Fixture
  : public ::testing::Test{

  using this_t	 = eigenDenseR9Fixture;
  using nat_mat_t = Eigen::MatrixXd;
  using mymat_t = pressio::containers::Matrix<nat_mat_t>;
  using nat_v_t	 = Eigen::VectorXd;
  using myvec_t	 = pressio::containers::Vector<nat_v_t>;

  // gold solution
  pressio::qr::test::qrGoldr9c4Sol<double> gold_;

  std::shared_ptr<mymat_t> A_;
  std::shared_ptr<myvec_t> v_;

  const int numRows_ = pressio::qr::test::numRows_;
  const int numVectors_ = pressio::qr::test::numVectors_;
  int numGlobalEntries_ = {};
  int localSize_ = {};

  virtual void SetUp(){
    numGlobalEntries_ = pressio::qr::test::numRows_;
    A_ = std::make_shared<mymat_t>(numRows_, numVectors_);
    v_ = std::make_shared<myvec_t>(numRows_);
    fillVector();
    fillMatrix();
  }//setup

  void fillMatrix(){
    (*A_)(0,0) = 3.2; (*A_)(0,1) = 1.2;  (*A_)(0,2) = 1.;
    (*A_)(1,0) = 1.2; (*A_)(1,2) = -2.2;
    (*A_)(2,1) = 4.0; (*A_)(2,3) = -2.;
    (*A_)(3,1) = 4.;
    (*A_)(4,2) = -1.; (*A_)(4,3) = -4.;
    (*A_)(5,0) = 0.2; (*A_)(5,1) = 5.;  (*A_)(5,2) = 1.;
    (*A_)(6,0) = 1.; (*A_)(6,1) = 1.1; (*A_)(6,2) = 1.25; (*A_)(6,3) = -3.;
    (*A_)(7,2) = 1.; (*A_)(7,1) = 0.1111; (*A_)(7,3) = 6.;
    std::cout << *(A_->data()) << std::endl;
  }

  void fillVector(){
    pressio::containers::ops::fill(*v_.get(), 1.);
  }

  void checkQFactor(const nat_mat_t & Q){
    EXPECT_EQ( Q.rows(), ::pressio::qr::test::numRows_);
    EXPECT_EQ( Q.cols(), ::pressio::qr::test::numVectors_);
    for (auto i=0; i<::pressio::qr::test::numRows_; i++){
      for (auto j=0; j<::pressio::qr::test::numVectors_; j++){
	EXPECT_NEAR( std::abs(Q(i,j)),
		     std::abs(gold_.trueQ_(i,j)), 1e-6);
      }
    }
  }

  virtual void TearDown(){}
};
// --------------------------------------------


#ifdef PRESSIO_ENABLE_TPL_TRILINOS
struct epetraR9Fixture
  : public ::testing::Test{

  using this_t	 = epetraR9Fixture;
  using nat_mv_t = Epetra_MultiVector;
  using mymvec_t = pressio::containers::MultiVector<nat_mv_t>;
  using nat_v_t	 = Epetra_Vector;
  using myvec_t	 = pressio::containers::Vector<nat_v_t>;

  // gold solution
  pressio::qr::test::qrGoldr9c4Sol<double> gold_;

  std::shared_ptr<Epetra_MpiComm> comm_;
  std::shared_ptr<Epetra_Map> rowMap_;
  std::shared_ptr<mymvec_t> A_;
  std::shared_ptr<myvec_t> v_;

  const int numVectors_ = pressio::qr::test::numVectors_;
  int rank_ = {};
  int numProc_ = {};
  int numGlobalEntries_ = {};
  int localSize_ = {};
  int shift_	 = {};

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = std::make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    rank_ = comm_->MyPID();
    numProc_ = comm_->NumProc();
    localSize_ = (rank_==0) ? 5 : 4;
    shift_ = (rank_==0) ? 0 : 5;
    assert(numProc_ == 2);
    numGlobalEntries_ = pressio::qr::test::numRows_;

    rowMap_ = std::make_shared<Epetra_Map>(numGlobalEntries_, 0, *comm_);
    A_ = std::make_shared<mymvec_t>(*rowMap_, numVectors_);
    v_ = std::make_shared<myvec_t>(*rowMap_);

    fillVector();
    fillMatrix();
  }

  void fillMatrix(){
    if(rank_==0){
      (*A_)(0,0) = 3.2; (*A_)(0,1) = 1.2;  (*A_)(0,2) = 1.;
      (*A_)(1,0) = 1.2; (*A_)(1,2) = -2.2;
      (*A_)(2,1) = 4.0; (*A_)(2,3) = -2.;
      (*A_)(3,1) = 4.;
      (*A_)(4,2) = -1.; (*A_)(4,3) = -4.;
    }
    if(rank_==1){
      (*A_)(0,0) = 0.2; (*A_)(0,1) = 5.;     (*A_)(0,2) = 1.;
      (*A_)(1,0) = 1.;  (*A_)(1,1) = 1.1;    (*A_)(1,2) = 1.25; (*A_)(1,3) = -3.;
      (*A_)(2,2) = 1.;  (*A_)(2,1) = 0.1111; (*A_)(2,3) = 6.;
    }
    //A_->data()->Print(std::cout);
  }

  void fillVector(){
    pressio::containers::ops::fill(*v_.get(), 1.);
  }

  void checkQFactor(const mymvec_t & Q){
    EXPECT_EQ( Q.extent(0), ::pressio::qr::test::numRows_);
    EXPECT_EQ( Q.numVectors(), ::pressio::qr::test::numVectors_);

    for (auto i=0; i<localSize_; i++)
      for (auto j=0; j<Q.numVectorsLocal(); j++){
	EXPECT_NEAR( std::abs(Q(i,j)),
		     std::abs(gold_.trueQ_(i+shift_,j)),
		     1e-6);
      }
  }

  virtual void TearDown(){}
};
// --------------------------------------------


struct tpetraR9Fixture
  : public ::testing::Test{

  using this_t	    = tpetraR9Fixture;
  using tcomm	    = Teuchos::Comm<int>;
  using map_t	    = Tpetra::Map<>;
  using nat_mvec_t  = Tpetra::MultiVector<>;
  using ST	    = typename nat_mvec_t::scalar_type;
  using LO	    = typename nat_mvec_t::local_ordinal_type;
  using GO	    = typename nat_mvec_t::global_ordinal_type;
  using nat_vec_t   = Tpetra::Vector<>;

  using mymvec_t    = pressio::containers::MultiVector<nat_mvec_t>;
  using mv_device_t = typename pressio::containers::details::traits<mymvec_t>::device_t;
  using myvec_t     = pressio::containers::Vector<nat_vec_t>;

  // gold solution
  pressio::qr::test::qrGoldr9c4Sol<ST> gold_;

  const int numVectors_ = pressio::qr::test::numVectors_;
  int rank_ = {};
  int numProc_ = {};
  int localSize_ = {};
  int shift_	 = {};

  int numGlobalEntries_;
  Teuchos::RCP<const tcomm> comm_;
  Teuchos::RCP<const map_t> contigMap_;
  std::shared_ptr<mymvec_t> A_;
  std::shared_ptr<myvec_t> v_;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    rank_ = comm_->getRank();
    numProc_ = comm_->getSize();
    localSize_ = (rank_==0) ? 5 : 4;
    shift_ = (rank_==0) ? 0 : 5;
    assert(numProc_==2);

    numGlobalEntries_ = pressio::qr::test::numRows_;
    contigMap_ = Teuchos::rcp(new map_t(numGlobalEntries_,0,comm_));
    A_ = std::make_shared<mymvec_t>(contigMap_, this_t::numVectors_);
    v_ = std::make_shared<myvec_t>(contigMap_);
    fillVector();
    fillMatrix();
  }

  void fillMatrix(){
    // get trilinos tpetra multivector object
    auto trilD = A_->data();
    trilD->sync<Kokkos::HostSpace>();

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
      v2d(0,0) = 3.2; v2d(0,1) = 1.2;  v2d(0,2) = 1.;
      v2d(1,0) = 1.2; v2d(1,2) = -2.2;
      v2d(2,1) = 4.0; v2d(2,3) = -2.;
      v2d(3,1) = 4.;
      v2d(4,2) = -1.; v2d(4,3) = -4.;
    }
    if(rank_==1){
      v2d(0,0) = 0.2; v2d(0,1) = 5.;     v2d(0,2) = 1.;
      v2d(1,0) = 1.;  v2d(1,1) = 1.1;    v2d(1,2) = 1.25; v2d(1,3) = -3.;
      v2d(2,2) = 1.;  v2d(2,1) = 0.1111; v2d(2,3) = 6.;
    }
    // sync from host to device
    trilD->sync<mv_device_t>();
  }

  void fillVector(){
    pressio::containers::ops::fill(*v_.get(), 1.);
  }

  void checkQFactor(const mymvec_t & Q)
  {
    EXPECT_EQ( Q.extent(0), ::pressio::qr::test::numRows_);
    EXPECT_EQ( Q.numVectors(), ::pressio::qr::test::numVectors_);
    for (auto j=0; j<Q.numVectorsLocal(); j++){
      auto colData = Q.data()->getData(j);
      for (auto i=0; i<localSize_; i++){
	EXPECT_NEAR( std::abs(colData[i]),
		     std::abs(gold_.trueQ_(i+shift_,j)),
		     1e-6);
      }
    }
  }

  virtual void TearDown(){}
};
#endif


#endif /* QR_FIXTURES_HPP_ */
