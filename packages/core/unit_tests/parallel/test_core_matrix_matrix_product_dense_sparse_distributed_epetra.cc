
#include <gtest/gtest.h>
#include "CORE_ALL"
#include "Epetra_MpiComm.h"

struct core_matrix_dense_sparse_distributed_epetraFix
  : public ::testing::Test{
public:
  int rank_;
  Epetra_MpiComm * Comm_;
  int MyPID_;
  int NumProc_;

  int nRowsDM_;
  int nColsDM_;
  Epetra_Map * dmMap_;
  core::Matrix<Epetra_CrsMatrix> * A_;

  int nRowsSM_;
  int nColsSM_;
  Epetra_Map * smMap_;
  core::Matrix<Epetra_MultiVector> * B_;
  
  virtual void SetUp()
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    Comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
    MyPID_ = Comm_->MyPID();
    NumProc_ = Comm_->NumProc();

    nRowsSM_ = NumProc_ * 3;
    smMap_ = new Epetra_Map(nRowsSM_, 0, *Comm_);
    A_ = new core::Matrix<Epetra_CrsMatrix>(*smMap_, 5);

    nRowsDM_ = A_->globalCols();
    nColsDM_ = 6;
    dmMap_ = new Epetra_Map(nRowsDM_, 0, *Comm_);
    B_ = new core::Matrix<Epetra_MultiVector>(*dmMap_, nColsDM_);
  }
  
  virtual void TearDown(){
    delete Comm_;
    delete dmMap_;
    delete smMap_;
    delete A_;
    delete B_;
  }
};

TEST_F(core_matrix_dense_sparse_distributed_epetraFix, Test1)
{

  //-----------
  // FILL A
  //-----------
  {
    int myNR = smMap_->NumMyElements();

    std::vector<int> mygid(myNR);
    smMap_->MyGlobalElements( mygid.data() );
    std::array<double,3> vals;
    std::array<int,3> colind;
    for (auto const & it : mygid){
      if(it == 0){
      	vals = {1., 2.};
      	colind = {0, 2};
      	A_->insertGlobalValues(it, 2, vals.data(), colind.data());
      }
      if(it == 2){
      	vals = {3., 4.};
      	colind = {0, 1};
      	A_->insertGlobalValues(it, 2, vals.data(), colind.data());
      }
    }
    A_->fillingIsCompleted();
    A_->data()->Print(std::cout);
  }

  //-----------
  // FILL B
  //-----------
  {
    int myNRows = dmMap_->NumMyElements();

    if (MyPID_ == 0){
      (*B_)(0,0) = 2.; (*B_)(0,1) = 1.;
      (*B_)(0,3) = 2.; (*B_)(0,4) = 2.;
      (*B_)(0,5) = 3.;

      (*B_)(2,0) = 3.;
    }
    B_->data()->Print(std::cout);
  }
  
  //-----------
  // product
  //-----------
  auto CC = core::matrixMatrixProduct(*A_, *B_);
  CC.data()->Print(std::cout);
  
}
