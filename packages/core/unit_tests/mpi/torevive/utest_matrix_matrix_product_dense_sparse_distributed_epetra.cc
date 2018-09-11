
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

  core::Matrix<Epetra_CrsMatrix> * A_;
  core::Matrix<Epetra_MultiVector> * B_;

  int nRowsSM_;
  int nColsSM_;
  Epetra_Map * smMap_;
  int nRowsDM_;
  int nColsDM_;
  Epetra_Map * dmMap_;
  
  virtual void SetUp()
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    Comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
    MyPID_ = Comm_->MyPID();
    NumProc_ = Comm_->NumProc();
    assert( NumProc_ == 3 );

    nRowsSM_ = 9;
    smMap_ = new Epetra_Map(nRowsSM_, 0, *Comm_);
    A_ = new core::Matrix<Epetra_CrsMatrix>(*smMap_, 5);

    nRowsDM_ = 13; //A_->globalCols();
    nColsDM_ = 5;
    dmMap_ = new Epetra_Map(nRowsDM_, 0, *Comm_);
    dmMap_->Print(std::cout);
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
    std::array<double,5> vals;
    std::array<int,5> colind;
    for (auto const & it : mygid){
      if(it == 0){
      	vals = {1., 1., 1., 3., 3.};
      	colind = {0, 4, 7, 11, 12};
      	A_->insertGlobalValues(it, 5, vals.data(), colind.data());
      }
      if(it == 4){
      	vals = {1., 3., 6.};
      	colind = {1, 3, 6};
      	A_->insertGlobalValues(it, 3, vals.data(), colind.data());
      }

    }
    A_->fillingIsCompleted(*dmMap_, *smMap_);
    A_->data()->Print(std::cout);
  }

  //-----------
  // FILL B
  //-----------
  {
    //    int myNRows = dmMap_->NumMyElements();
    if (MyPID_ == 0){
      (*B_)(0,0) = 1.; (*B_)(0,2) = 1.;
      (*B_)(0,3) = 2.; (*B_)(0,4) = 2.;

      (*B_)(1,2) = 2.;
      (*B_)(3,2) = 4.;
      (*B_)(4,0) = 3.;
    }

    if (MyPID_ == 1){
      (*B_)(2,0) = 1.;
      (*B_)(1,2) = 2.;
    }

    if (MyPID_ == 2){
      (*B_)(2,0) = 2.;
      (*B_)(3,0) = 2.;
    }

    B_->data()->Print(std::cout);
  }
  
  // //------------------
  // // product: C = A B
  // //------------------
  // auto CC = core::mat_ops::product(*A_, *B_);
  // CC.data()->Print(std::cout);

  // assert( CC.globalRows() == 9);
  // assert( CC.globalCols() == 5);
  // static_assert( std::is_same<decltype(CC), core::Matrix<Epetra_MultiVector>>::value, "" );
  // if (MyPID_ == 0){
  //   EXPECT_DOUBLE_EQ( CC(0,0), 17. );
  //   EXPECT_DOUBLE_EQ( CC(0,1), 0. );
  //   EXPECT_DOUBLE_EQ( CC(0,2), 1. );
  //   EXPECT_DOUBLE_EQ( CC(0,3), 2. );
  //   EXPECT_DOUBLE_EQ( CC(0,4), 2. );
  // }
  // if (MyPID_ == 1)
  //   EXPECT_DOUBLE_EQ( CC(1,2), 26. );
  
}
