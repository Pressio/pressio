
#include <gtest/gtest.h>
#include "Epetra_MpiComm.h"
#include "CONTAINERS_VECTOR"
#include "CONTAINERS_MULTI_VECTOR"
#include "CONTAINERS_MATRIX"

struct containers_matrix_dense_distributed_epetraFix
  : public ::testing::Test{
public:
  int rank_;
  Epetra_MpiComm * Comm_;
  int MyPID_;
  int NumProc_;
  const int localSize_ = 5;
  const int numVectors_ = 6;
  int numGlobalEntries_;
  Epetra_Map * dataMap_;
  Epetra_MultiVector * x_;
  
  virtual void SetUp()
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    Comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
    MyPID_ = Comm_->MyPID();
    NumProc_ = Comm_->NumProc();
    numGlobalEntries_ = Comm_->NumProc() * localSize_;
    dataMap_ = new Epetra_Map(numGlobalEntries_, 0, *Comm_);
    x_ = new Epetra_MultiVector(*dataMap_, numVectors_);
  }

  Epetra_MultiVector & getNative(){
    return *x_;
  }

  int getRank() const{
    return rank_;
  }
  int getNumProc() const{
    return NumProc_;
  }
  int numGlobalEntries() const{
    return numGlobalEntries_;
  }
  int numLocalEntries() const{
    return localSize_;
  }
  const Epetra_Map * getMap() const{
    return dataMap_;
  }    
  virtual void TearDown(){
    delete Comm_;
    delete dataMap_;
    delete x_;
  }
};

TEST_F(containers_matrix_dense_distributed_epetraFix, EpetraDenseDistConstructor)
{
  using nat_t = Epetra_MultiVector;
  using my_t = containers::Matrix<nat_t>;
  STATIC_ASSERT_IS_CONTAINERS_MATRIX_WRAPPER(my_t);

  my_t b( *getMap(), numVectors_ );
  ASSERT_FALSE( b.empty() );

  EXPECT_EQ( b.globalCols(), 6 );
  EXPECT_EQ( b.localCols(), 6 );
  EXPECT_EQ( b.globalRows(), 15 );
  EXPECT_EQ( b.localRows(), 5);

  for (int i=0; i<b.localRows(); i++)
    for (int j=0; j<b.globalCols(); j++)
      EXPECT_NEAR( 0.0, b(i,j), 1e-12);

  if(MyPID_==0)
    b.replaceGlobalValue(1,1, 43.3);
  if(MyPID_==0)
    EXPECT_NEAR( 43.3, b(1,1), 1e-12);
  else
    EXPECT_NEAR( 0.0, b(1,1), 1e-12);

  b.setZero();
  for (int i=0; i<b.localRows(); i++)
    for (int j=0; j<b.globalCols(); j++)
      EXPECT_NEAR( 0.0, b(i,j), 1e-12);  

}


// TEST_F(containers_matrix_dense_distributed_epetraFix, transpose)
// {
//   using nat_t = Epetra_MultiVector;
//   using my_t = containers::Matrix<nat_t>;
//   STATIC_ASSERT_IS_CONTAINERS_MATRIX_WRAPPER(my_t);

//   my_t A( *getMap(), numVectors_ );
//   auto CC = containers::mat_ops::transpose(A);
//   assert( CC.globalRows() == A.globalCols() );
//   assert( CC.globalCols() == A.globalRows() );

//   // missing test here because we are missing implementation
  
// }
