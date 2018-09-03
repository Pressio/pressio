
#include <gtest/gtest.h>
#include "vector/core_vector_meta.hpp"
#include "vector/core_vector_distributed_epetra.hpp"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"
#include "Epetra_MpiComm.h"


// struct core_vector_distributed_epetraFix
//   : public ::testing::Test{
// public:

//   int rank_;
//   Epetra_MpiComm * Comm_;
//   int MyPID_;
//   int NumProc_;
//   const int localSize_ = 5;
//   int numGlobalEntries_;
//   Epetra_Map * contigMap_;
//   Epetra_Vector * x_;
  
//   virtual void SetUp()
//   {
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
//     Comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
//     MyPID_ = Comm_->MyPID();
//     NumProc_ = Comm_->NumProc();
//     numGlobalEntries_ = Comm_->NumProc() * localSize_;
//     contigMap_ = new Epetra_Map(numGlobalEntries_, 0, *Comm_);
//     x_ = new Epetra_Vector(*contigMap_);
//   }

//   int getRank() const{
//     return rank_;
//   }
//   int getNumProc() const{
//     return NumProc_;
//   }
//   int numGlobalEntries() const{
//     return numGlobalEntries_;
//   }
//   int numLocalEntries() const{
//     return localSize_;
//   }
//   const Epetra_Map * getMap() const{
//     return contigMap_;
//   }    
//   const Epetra_Vector * getVector() const{
//     return x_;
//   }
//   void fillScalar(double value){
//     x_->PutScalar(value);
//   }
    
//   virtual void TearDown(){
//     delete Comm_;
//     delete contigMap_;
//     delete x_;
//   }
// };



// using core_vector_distributed_epetra_DeathTest
// = core_vector_distributed_epetraFix;

// TEST_F(core_vector_distributed_epetra_DeathTest,
//        EpetraVectorSubscriptOperator)
// {
//   if (getRank()==0){
//     using myvec_t = core::Vector<Epetra_Vector>;
//     myvec_t v1( *getMap() );
//     int localSize = numLocalEntries();
//     //ASSERT_DEATH(v1[localSize+1]=4.0, "Assertion failed:");
//     //...
//   }
// }
