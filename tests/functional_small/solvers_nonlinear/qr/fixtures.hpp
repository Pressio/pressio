
#ifndef QR_UTEST_FIXTURES_HPP_
#define QR_UTEST_FIXTURES_HPP_

#include <gtest/gtest.h>
#include "pressio/solvers_nonlinear.hpp"
#include "qr_r9c4_gold.hpp"

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#ifdef PRESSIO_ENABLE_EPETRA
#include "Epetra_MpiComm.h"
#endif // PRESSIO_ENABLE_EPETRA
#include "Eigen/Dense"
#endif // PRESSIO_ENABLE_TPL_TRILINOS

namespace pressio{ namespace qr{ namespace test{
  static constexpr int numVectors_ = 4;
  static constexpr int numRows_ = 9;
}}}// end namespace pressio::qr::test


struct eigenDenseR9Fixture
  : public ::testing::Test{

  using this_t	 = eigenDenseR9Fixture;
  using matrix_type = Eigen::MatrixXd;
  using vector_type  = Eigen::VectorXd;
  pressio::qr::test::qrGoldr9c4Sol<double> gold_;
  matrix_type A_;
  vector_type v_;

  const int numRows_ = pressio::qr::test::numRows_;
  const int numVectors_ = pressio::qr::test::numVectors_;
  int numGlobalEntries_ = {};
  int localSize_ = {};

  virtual void SetUp(){
    numGlobalEntries_ = pressio::qr::test::numRows_;
    A_.resize(numRows_, numVectors_);
    v_.resize(numRows_);
    A_.setZero();
    v_.setZero();
    fillVector();
    fillMatrix();
  }//setup

  void fillMatrix()
  {
    A_(0,0) = 3.2; A_(0,1) = 1.2;  A_(0,2) = 1.;
    A_(1,0) = 1.2; A_(1,2) = -2.2;
    A_(2,1) = 4.0; A_(2,3) = -2.;
    A_(3,1) = 4.;
    A_(4,2) = -1.; A_(4,3) = -4.;
    A_(5,0) = 0.2; A_(5,1) = 5.;  A_(5,2) = 1.;
    A_(6,0) = 1.; A_(6,1) = 1.1; A_(6,2) = 1.25; A_(6,3) = -3.;
    A_(7,2) = 1.; A_(7,1) = 0.1111; A_(7,3) = 6.;
    // std::cout << A_ << std::endl;
  }

  void fillVector(){
    pressio::ops::fill(v_, 1.);
  }

  template<class T>
  void checkQFactor(const T & Q)
  {
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
#ifdef PRESSIO_ENABLE_EPETRA
struct epetraR9Fixture
  : public ::testing::Test{

  using this_t	 = epetraR9Fixture;
  using mymvec_t = Epetra_MultiVector;
  using myvec_t	 = Epetra_Vector;

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
    A_->PutScalar(0);
    v_->PutScalar(0);

    fillVector();
    fillMatrix();
  }

  void fillMatrix()
  {
    if(rank_==0){
      (*A_)[0][0] = 3.2; (*A_)[1][0] = 1.2;  (*A_)[2][0] = 1.;
      (*A_)[0][1] = 1.2; (*A_)[2][1] = -2.2;
      (*A_)[1][2] = 4.0; (*A_)[3][2] = -2.;
      (*A_)[1][3] = 4.;
      (*A_)[2][4] = -1.; (*A_)[3][4] = -4.;
    }
    if(rank_==1){
      (*A_)[0][0] = 0.2; (*A_)[1][0] = 5.;     (*A_)[2][0] = 1.;
      (*A_)[0][1] = 1.;  (*A_)[1][1] = 1.1;    (*A_)[2][1] = 1.25; (*A_)[3][1] = -3.;
      (*A_)[2][2] = 1.;  (*A_)[1][2] = 0.1111; (*A_)[3][2] = 6.;
    }
    // A_->Print(std::cout);
  }

  void fillVector(){
    pressio::ops::fill(*v_.get(), 1.);
  }

  void checkQFactor(const mymvec_t & Q){
    EXPECT_EQ( Q.GlobalLength(),   ::pressio::qr::test::numRows_);
    EXPECT_EQ( Q.NumVectors(), ::pressio::qr::test::numVectors_);

    for (auto i=0; i<localSize_; i++)
      for (auto j=0; j<Q.NumVectors(); j++){
	       EXPECT_NEAR(std::abs(Q[j][i]), std::abs(gold_.trueQ_(i+shift_,j)), 1e-6);
      }
  }

  virtual void TearDown(){}
};
#endif // PRESSIO_ENABLE_EPETRA
// --------------------------------------------
struct tpetraR9Fixture
  : public ::testing::Test{

  using this_t	    = tpetraR9Fixture;
  using tcomm	    = Teuchos::Comm<int>;
  using map_t	    = Tpetra::Map<>;

  using mymvec_t    = Tpetra::MultiVector<>;
  using myvec_t     = Tpetra::Vector<>;
  using ST      = typename mymvec_t::scalar_type;
  using LO      = typename mymvec_t::local_ordinal_type;
  using GO      = typename mymvec_t::global_ordinal_type;

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
    A_->putScalar(0);
    v_->putScalar(0);

    fillVector();
    fillMatrix();
  }

  void fillMatrix()
  {
    // get trilinos tpetra multivector object
    auto trilD = A_;
    // trilD->sync<Kokkos::HostSpace>();

    /*--------------------------------------------
     * (1): modify the host view and then sync
     * most likely, host and device will be same unless we run CUDA
     * so in theory we should not worry about syncing but it
     * does not hurt to do it anyway
     //--------------------------------------------*/
    auto v2d = trilD->getLocalView<Kokkos::HostSpace>(Tpetra::Access::ReadWriteStruct());
    auto c0 = Kokkos::subview(v2d, Kokkos::ALL(), 0);
    // //we are going to change the host view
    // trilD->modify<Kokkos::HostSpace>();

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
    // // sync from host to device
    // trilD->sync<mv_device_t>();
  }

  void fillVector(){
    pressio::ops::fill(*v_.get(), 1.);
  }

  void checkQFactor(const mymvec_t & Q)
  {
    EXPECT_EQ( Q.getGlobalLength(), ::pressio::qr::test::numRows_);
    EXPECT_EQ( Q.getNumVectors(), ::pressio::qr::test::numVectors_);
    for (int j=0; j< (int)Q.getNumVectors(); j++){
      auto colData = Q.getData(j);
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
