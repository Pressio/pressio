
#ifndef CONTAINERS_FIXTURES_EPETRA_ONLY_FIXTURES_HPP_
#define CONTAINERS_FIXTURES_EPETRA_ONLY_FIXTURES_HPP_

#include <gtest/gtest.h>
#include "Epetra_MpiComm.h"
#include "pressio_containers.hpp"

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
