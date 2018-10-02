
#ifndef CORE_FIXTURES_EPETRA_ONLY_FIXTURES_HPP_
#define CORE_FIXTURES_EPETRA_ONLY_FIXTURES_HPP_

#include <gtest/gtest.h>
#include "Epetra_MpiComm.h"
#include "CORE_ALL"

struct epetraVectorGlobSize15Fixture
  : public ::testing::Test{

public:
  Epetra_MpiComm * comm_;
  int rank_;
  int numProc_;
  const int localSize_ = 5;
  int numGlobalEntries_;
  Epetra_Map * contigMap_;
  Epetra_Vector * x_;
  
  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
    rank_ = comm_->MyPID();
    numProc_ = comm_->NumProc();
    assert(numProc_==3);
    
    numGlobalEntries_ = numProc_ * localSize_;
    contigMap_ = new Epetra_Map(numGlobalEntries_, 0, *comm_);
    x_ = new Epetra_Vector(*contigMap_);
  }
  
  virtual void TearDown(){
    delete comm_;
    delete contigMap_;
    delete x_;
  }
};
//-----------------------------------------------------------


struct epetraMultiVectorR9C4VecS9Fixture
  : public ::testing::Test{
  
public:
  int rank_;
  Epetra_MpiComm * comm_;
  int numProc_;
  const int localSize_ = 3;
  const int numVectors_ = 4;
  int numGlobalEntries_;
  Epetra_Map * dataMap_;
  Epetra_MultiVector * mv_;
  Epetra_Vector * x_;
  
  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
    rank_ = comm_->MyPID();
    numProc_ = comm_->NumProc();
    assert(numProc_ == 3);
    
    numGlobalEntries_ = numProc_ * localSize_;
    dataMap_ = new Epetra_Map(numGlobalEntries_, 0, *comm_);
    mv_ = new Epetra_MultiVector(*dataMap_, numVectors_);
    x_ = new Epetra_Vector(*dataMap_);
  }

  virtual void TearDown(){
    delete comm_;
    delete dataMap_;
    delete mv_;
    delete x_;
  }
};
//-----------------------------------------------------------


#endif /* CORE_FIXTURES_EPETRA_ONLY_FIXTURES_HPP_ */

