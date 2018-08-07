
#include <gtest/gtest.h>
#include "CORE_ALL"
#include "Epetra_LocalMap.h"
#include "Epetra_MpiComm.h"

struct core_matrix_dense_distributed_epetraFix
  : public ::testing::Test{
public:
  int rank_;
  Epetra_MpiComm * Comm_;
  int MyPID_;
  int NumProc_;
  Epetra_Map * contigMapA_;
  Epetra_Map * contigMapB_;
  core::Matrix<Epetra_MultiVector> * A_;
  core::Matrix<Epetra_MultiVector> * B_;
  
  virtual void SetUp()
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    Comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
    MyPID_ = Comm_->MyPID();
    NumProc_ = Comm_->NumProc();

    assert( NumProc_ == 3);
    
    // A
    int nRowsA = 3 * NumProc_;
    int nRowsB = 2 * NumProc_;
    int nColsA = nRowsB;
    int nColsB = 4;
        
    contigMapA_ = new Epetra_Map(nRowsA, 0, *Comm_);
    contigMapB_ = new Epetra_Map(nRowsB, 0, *Comm_);
    
    A_ = new core::Matrix<Epetra_MultiVector>(*contigMapA_, nColsA);
    B_ = new core::Matrix<Epetra_MultiVector>(*contigMapB_, nColsB);
  }
  
  virtual void TearDown(){
    delete Comm_;
    delete contigMapA_;
    delete contigMapB_;
    delete A_;
    delete B_;
  }
};


TEST_F(core_matrix_dense_distributed_epetraFix, Test1)
{
  //-----------
  // FILL A
  // 9 x 6 matrix 
  //-----------
  {
    if (MyPID_ == 0){
      (*A_)(0,0) = 2.; (*A_)(0,1) = 1.;
      (*A_)(0,3) = 2.; (*A_)(0,4) = 2.;
      (*A_)(0,5) = 3.;

      (*A_)(2,1) = 2.; (*A_)(2,3) = 3.;
    }
    if (MyPID_ == 1){
      (*A_)(2,0) = 1.; (*A_)(2,2) = 1.;
      (*A_)(2,4) = 1.; (*A_)(2,5) = 1.;
    }
    if (MyPID_ == 2){
      (*A_)(1,1) = 2.; (*A_)(1,3) = 2.;
      (*A_)(1,4) = 2.; (*A_)(1,5) = 2.;
    }

    A_->data()->Print(std::cout);
  }

  //-----------
  // FILL B
  // 6 x 4 matrix 
  //-----------
  {
    if (MyPID_ == 0){
      (*B_)(0,0) = 2.; (*B_)(1,0) = 1.;
      (*B_)(1,2) = 1.; (*B_)(0,3) = 2.;
    }
    if (MyPID_ == 1){
      (*B_)(0,0) = 2.; (*B_)(0,3) = 1.;
      (*B_)(1,0) = 2.; (*B_)(1,2) = 3.;
    }
    if (MyPID_ == 2){
      (*B_)(0,3) = 3.;
    }
    B_->data()->Print(std::cout);
  }

  //-----------
  // product C
  // 9 x 4 matrix 
  //-----------
  auto CC = core::matrixMatrixProduct(*A_, *B_);
  assert( CC.globalRows() == 9 );
  assert( CC.globalCols() == 4 );
  CC.data()->Print(std::cout);

  if (MyPID_ == 0){
    EXPECT_DOUBLE_EQ( CC(0,0), 9.0 );
    EXPECT_DOUBLE_EQ( CC(0,2), 7.0 );
    EXPECT_DOUBLE_EQ( CC(0,3), 10.0 );
    EXPECT_DOUBLE_EQ( CC(2,0), 8.0 );
    EXPECT_DOUBLE_EQ( CC(2,2), 11.0 );
  }
  if (MyPID_ == 1){
    EXPECT_DOUBLE_EQ( CC(2,0), 4.0 );
    EXPECT_DOUBLE_EQ( CC(2,3), 6.0 );
  }
  if (MyPID_ == 2){
    EXPECT_DOUBLE_EQ( CC(1,0), 6.0 );
    EXPECT_DOUBLE_EQ( CC(1,2), 8.0 );
  }  

}
