
#include <gtest/gtest.h>
#include "CORE_ALL"
#include "Epetra_MpiComm.h"

struct core_matrix_vec_product_distributed_epetraFix
  : public ::testing::Test{
public:
  int rank_;
  Epetra_MpiComm * Comm_;
  int MyPID_;
  int NumProc_;

  core::Matrix<Epetra_MultiVector> * Ad_;
  core::Vector<Epetra_Vector> * b2_;
  int nRowsDM_;
  Epetra_Map * dmMap_;
  Epetra_Map * b2Map_;
  
  virtual void SetUp()
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    Comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
    MyPID_ = Comm_->MyPID();
    NumProc_ = Comm_->NumProc();
    assert( NumProc_ == 3 );

    nRowsDM_ = 9;
    dmMap_ = new Epetra_Map(nRowsDM_, 0, *Comm_);
    Ad_ = new core::Matrix<Epetra_MultiVector>(*dmMap_, 5);
    b2Map_ = new Epetra_Map(5, 0, *Comm_);
    b2Map_->Print(std::cout);
    b2_ = new core::Vector<Epetra_Vector>(*b2Map_);
  }
  
  virtual void TearDown(){
    delete Comm_;
    delete b2Map_;
    delete dmMap_;
    delete Ad_;
    delete b2_;
  }
};


TEST_F(core_matrix_vec_product_distributed_epetraFix, DenseMatTimesVector){

  //-----------
  // FILL A
  //-----------
  {
    if (MyPID_ == 0){
      (*Ad_)(0,0) = 1.0;
      (*Ad_)(0,2) = 2.0;
      (*Ad_)(0,3) = 3.0;
      (*Ad_)(0,4) = 1.0;
    }
    if (MyPID_ == 1){
      (*Ad_)(1,0) = 3.0;
      (*Ad_)(1,2) = 2.0;
      (*Ad_)(1,3) = 3.0;
      (*Ad_)(1,4) = 4.0;
    }

    Ad_->data()->Print(std::cout);
  }

  //-----------
  // FILL b
  //-----------
  {
    if (MyPID_ == 0){
      (*b2_)[0] = 1.;
      (*b2_)[1] = 2.;
    }
    if (MyPID_ == 1){
      (*b2_)[0] = 1.;
    }
    if (MyPID_ == 2){
      (*b2_)[0] = 1.;
    }
    b2_->data()->Print(std::cout);
  }
  
  //------------------
  // product: b = A b
  //------------------
  auto c = core::mat_ops::product( *Ad_, *b2_ );  
  c.data()->Print(std::cout);

  assert( c.globalSize() == 9);
  static_assert( std::is_same<decltype(c),
  		 core::Vector<Epetra_Vector>>::value, "" );
  if (MyPID_ == 0){
    EXPECT_DOUBLE_EQ( c[0], 4. );
  }
  if (MyPID_ == 1)
    EXPECT_DOUBLE_EQ( c[1], 9. );
  
}
