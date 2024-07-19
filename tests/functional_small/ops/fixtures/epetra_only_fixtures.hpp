
#ifndef CONTAINERS_FIXTURES_EPETRA_ONLY_FIXTURES_HPP_
#define CONTAINERS_FIXTURES_EPETRA_ONLY_FIXTURES_HPP_

#include <gtest/gtest.h>
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"

struct epetraVectorGlobSize15Fixture
  : public ::testing::Test{

public:
  std::shared_ptr<Epetra_MpiComm> comm_;
  int rank_;
  int numProc_;
  const int localSize_ = 5;
  int numGlobalEntries_;
  std::shared_ptr<Epetra_Map> contigMap_;
  std::shared_ptr<Epetra_Vector> myVector_;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = std::make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    rank_ = comm_->MyPID();
    numProc_ = comm_->NumProc();
    EXPECT_EQ(numProc_,3);

    numGlobalEntries_ = numProc_ * localSize_;
    contigMap_ = std::make_shared<Epetra_Map>(numGlobalEntries_, 0, *comm_);
    myVector_ = std::make_shared<Epetra_Vector>(*contigMap_);
  }

  virtual void TearDown(){}
};
//-----------------------------------------------------------


struct epetraMultiVectorGlobSize15Fixture
  : public ::testing::Test{

public:
  std::shared_ptr<Epetra_MpiComm> comm_;
  int rank_;
  int numProc_;
  const int numVecs_ = 4;
  const int localSize_ = 5;
  int numGlobalEntries_;
  std::shared_ptr<Epetra_Map> contigMap_;
  std::shared_ptr<Epetra_Map> map_to_all_;
  std::shared_ptr<Epetra_Import> importer_;
  std::shared_ptr<Epetra_MultiVector> myMv_;
  std::shared_ptr<Epetra_Vector> x_epetra;
  std::shared_ptr<Epetra_Vector> y_epetra;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = std::make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    rank_ = comm_->MyPID();
    numProc_ = comm_->NumProc();
    EXPECT_EQ(numProc_,3);

    numGlobalEntries_ = numProc_ * localSize_;
    contigMap_ = std::make_shared<Epetra_Map>(numGlobalEntries_, 0, *comm_);
    // Note: this importer sends whole object to all ranks
    map_to_all_ = std::make_shared<Epetra_Map>(numGlobalEntries_, numGlobalEntries_, 0, *comm_);
    importer_ = std::make_shared<Epetra_Import>(*map_to_all_, *contigMap_);

    myMv_ = std::make_shared<Epetra_MultiVector>(*contigMap_, numVecs_);
    myMv_->PutScalar(1.);
    for (int i = 0; i < localSize_; ++i) {
      for (int j = 0; j < numVecs_; ++j) {
        // generate rank-unique int values
        (*myMv_)[j][i] = (double)((rank_ * localSize_ + i) * numVecs_ + j + 1.);
      }
    }
    x_epetra = std::make_shared<Epetra_Vector>(*contigMap_);
    for (int j = 0; j < localSize_; ++j) {
      // generate rank-unique int values
      (*x_epetra)[j] = (double)(rank_ * localSize_ + j + 1.);
    }
    y_epetra = std::make_shared<Epetra_Vector>(*contigMap_);
    y_epetra->PutScalar(3.);
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
  std::shared_ptr<Epetra_MultiVector> myMv_;
  std::shared_ptr<Epetra_Vector> myVector_;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = std::make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    rank_ = comm_->MyPID();
    numProc_ = comm_->NumProc();
    EXPECT_EQ(numProc_,3);

    numGlobalEntries_ = numProc_ * localSize_;
    dataMap_ = std::make_shared<Epetra_Map>(numGlobalEntries_, 0, *comm_);
    myMv_ = std::make_shared<Epetra_MultiVector>(*dataMap_, numVectors_);
    myVector_ = std::make_shared<Epetra_Vector>(*dataMap_);
  }

  virtual void TearDown(){}
};

#endif /* CONTAINERS_FIXTURES_EPETRA_ONLY_FIXTURES_HPP_ */
