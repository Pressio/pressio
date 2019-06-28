
#ifndef ALGEBRA_FIXTURES_EPETRA_ONLY_FIXTURES_HPP_
#define ALGEBRA_FIXTURES_EPETRA_ONLY_FIXTURES_HPP_

#include <gtest/gtest.h>
#include "Epetra_MpiComm.h"
#include "ALGEBRA_ALL"


struct epetraVectorGlobSize15Fixture
  : public ::testing::Test{

public:
  std::shared_ptr<Epetra_MpiComm> comm_;
  int rank_;
  int numProc_;
  const int localSize_ = 5;
  int numGlobalEntries_;
  std::shared_ptr<Epetra_Map> contigMap_;
  std::shared_ptr<Epetra_Vector> x_;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = std::make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    rank_ = comm_->MyPID();
    numProc_ = comm_->NumProc();
    EXPECT_EQ(numProc_,3);

    numGlobalEntries_ = numProc_ * localSize_;
    contigMap_ = std::make_shared<Epetra_Map>(numGlobalEntries_, 0, *comm_);
    x_ = std::make_shared<Epetra_Vector>(*contigMap_);
  }

  virtual void TearDown(){}
};
//-----------------------------------------------------------


struct epetraMultiVectorR9C4VecS9Fixture
  : public ::testing::Test{

public:
  int rank_;
  std::shared_ptr<Epetra_MpiComm> comm_;
  int numProc_;
  const int localSize_ = 3;
  const int numVectors_ = 4;
  int numGlobalEntries_;
  std::shared_ptr<Epetra_Map> dataMap_;
  std::shared_ptr<Epetra_MultiVector> mv_;
  std::shared_ptr<Epetra_Vector> x_;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = std::make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    rank_ = comm_->MyPID();
    numProc_ = comm_->NumProc();
    EXPECT_EQ(numProc_,3);

    numGlobalEntries_ = numProc_ * localSize_;
    dataMap_ = std::make_shared<Epetra_Map>(numGlobalEntries_, 0, *comm_);
    mv_ = std::make_shared<Epetra_MultiVector>(*dataMap_, numVectors_);
    x_ = std::make_shared<Epetra_Vector>(*dataMap_);
  }

  virtual void TearDown(){}
};
//-----------------------------------------------------------


struct epetraSparseMatR7MultiVectorR9C4Fixture
  : public ::testing::Test{

public:
  int rank_;
  std::shared_ptr<Epetra_MpiComm> comm_;
  int numProc_;

  // for crs matrix
  using mymat_w_t = rompp::algebra::Matrix<Epetra_CrsMatrix>;
  std::shared_ptr<Epetra_Map> dataMapSM_;
  int nRowsSM_ = 7;
  std::shared_ptr<mymat_w_t> sm_;

  // multivector
  using mymv_w_t =
    rompp::algebra::MultiVector<Epetra_MultiVector>;
  const int numVectors_ = 4;
  int nRowsMV_ = 9;
  std::shared_ptr<Epetra_Map> dataMapMV_;
  std::shared_ptr<mymv_w_t> mv_;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = std::make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    rank_ = comm_->MyPID();
    numProc_ = comm_->NumProc();
    EXPECT_EQ(numProc_,3);

    dataMapSM_ = std::make_shared<Epetra_Map>(nRowsSM_, 0, *comm_);
    sm_ = std::make_shared<mymat_w_t>(*dataMapSM_, 3);
    dataMapMV_ = std::make_shared<Epetra_Map>(nRowsMV_, 0, *comm_);
    mv_ = std::make_shared<mymv_w_t>(*dataMapMV_, numVectors_);
  }

  void fillMultiVector(){
    auto & mv = *mv_;
    if (rank_==0){
      mv(0,0)=0; mv(0,1)=1.; mv(0,2)=1.; mv(0,3)=0;
      mv(1,0)=0; mv(1,1)=0.; mv(1,2)=0.; mv(1,3)=1.;
      mv(2,0)=0; mv(2,1)=1.; mv(2,2)=0.; mv(2,3)=0;
    }

    if (rank_==1){
      mv(0,0)=1; mv(0,1)=1.; mv(0,2)=0.; mv(0,3)=0;
      mv(1,0)=0; mv(1,1)=2.; mv(1,2)=1.; mv(1,3)=0.;
      mv(2,0)=0; mv(2,1)=1.; mv(2,2)=0.; mv(2,3)=0;
    }

    if (rank_==2){
      mv(0,0)=2; mv(0,1)=2.; mv(0,2)=1.; mv(0,3)=0;
      mv(1,0)=0; mv(1,1)=1.; mv(1,2)=3.; mv(1,3)=0.;
      mv(2,0)=2; mv(2,1)=1.; mv(2,2)=1.; mv(2,3)=0;
    }
  }//end fillMultiVector


  void fillCrsMatrix(){

    // 1 0 2 0 0 3 0
    // 1 0 1 0 1 1 0
    // 0 0 1 2 1 0 0
    // 1 0 0 0 2 1 3
    // 0 0 1 1 0 0 0
    // 0 0 3 3 0 1 0
    // 0 0 0 3 0 1 0

    std::array<double,4> va;
    std::array<int,4> ci;
    if (rank_==0){
      //row 0
      va = {{1., 2., 3.}};  ci = {{0, 2, 5}};
      sm_->insertGlobalValues(0, 3, va.data(), ci.data());
      //row 1
      va = {{1.,1.,1.,1.}};  ci = {{0,2,4,5}};
      sm_->insertGlobalValues(1, 4, va.data(), ci.data());
      //row 2
      va = {{1.,2.,1.}};  ci = {{2,3,4}};
      sm_->insertGlobalValues(2, 3, va.data(), ci.data());
    }

    if (rank_==1){
      //row 3
      va = {{1.,2.,1,3.}};  ci = {{0,4,5,6}};
      sm_->insertGlobalValues(3, 4, va.data(), ci.data());
      //row 4
      va = {{1.,1.}};  ci = {{2,3}};
      sm_->insertGlobalValues(4, 2, va.data(), ci.data());
    }

    if (rank_==2){
      //row 5
      va = {{3.,3.,1}};  ci = {{2,3,5}};
      sm_->insertGlobalValues(5, 3, va.data(), ci.data());
      //row 6
      va = {{3,1.}};  ci = {{3,5}};
      sm_->insertGlobalValues(6, 2, va.data(), ci.data());
    }
  }//end fill

  virtual void TearDown(){}
};
//-----------------------------------------------------------


#endif /* ALGEBRA_FIXTURES_EPETRA_ONLY_FIXTURES_HPP_ */
