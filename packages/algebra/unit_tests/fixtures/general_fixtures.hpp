
#ifndef ALGEBRA_FIXTURES_GENERAL_FIXTURES_HPP_
#define ALGEBRA_FIXTURES_GENERAL_FIXTURES_HPP_

#include <gtest/gtest.h>
#include "ALGEBRA_ALL"

// struct algebra_general_fixture
//   : public ::testing::Test{
// public:
//   using namespace rompp;

//   int rank_;
//   Epetra_MpiComm * Comm_;
//   int MyPID_;
//   int NumProc_;
//   const int localSize_ = 5;
//   const int numVectors_ = 3;
//   int numGlobalEntries_;
//   Epetra_Map * dataMap_;
//   Epetra_MultiVector * x_;
  
//   virtual void SetUp()
//   {
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
//     Comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
//     MyPID_ = Comm_->MyPID();
//     NumProc_ = Comm_->NumProc();
//     numGlobalEntries_ = Comm_->NumProc() * localSize_;
//     dataMap_ = new Epetra_Map(numGlobalEntries_, 0, *Comm_);
//     x_ = new Epetra_MultiVector(*dataMap_, numVectors_);
//   }
//   virtual void TearDown(){
//     delete Comm_;
//     delete dataMap_;
//     delete x_;
//   }

// };//end struct
// //-------------------------------------------------------


#endif /* ALGEBRA_FIXTURES_GENERAL_FIXTURES_HPP_ */

