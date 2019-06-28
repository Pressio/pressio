
#include <gtest/gtest.h>
#include "Epetra_MpiComm.h"
#include "ALGEBRA_VECTOR"
#include "ALGEBRA_MATRIX"

struct algebra_matrix_sparse_distributed_epetraFix
  : public ::testing::Test{
public:
  int rank_;
  Epetra_MpiComm * Comm_;
  int MyPID_;
  int NumProc_;
  const int localSize_ = 5;
  int numGlobalEntries_;
  Epetra_Map * contigMap_;
  algebra::Matrix<Epetra_CrsMatrix> * A_;
  
  virtual void SetUp()
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    Comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
    MyPID_ = Comm_->MyPID();
    NumProc_ = Comm_->NumProc();
    numGlobalEntries_ = Comm_->NumProc() * localSize_;
    contigMap_ = new Epetra_Map(numGlobalEntries_, 0, *Comm_);
    A_ = new algebra::Matrix<Epetra_CrsMatrix>(*contigMap_, 3);
  }
  
  virtual void TearDown(){
    delete Comm_;
    delete contigMap_;
    delete A_;
  }
};


TEST_F(algebra_matrix_sparse_distributed_epetraFix,
       EpetraCrsMatrixTranspose)
{
  int myN = contigMap_->NumMyElements();

  // get my global indices
  std::vector<int> mygid(myN);
  contigMap_->MyGlobalElements( mygid.data() );

  std::array<double,3> vals;
  std::array<int,3> colind;
  
  for (auto const & it : mygid){
    if(it == 0){
      vals = {1., 5., 2.};
      colind = {0, 3, 4};
      A_->insertGlobalValues(it, 3, vals.data(), colind.data());
    }
    if(it == 3){
      vals = {11., 55., 22.};
      colind = {6, 8, 11};
      A_->insertGlobalValues(it, 3, vals.data(), colind.data());
    }
    if(it == 11){
      vals = {111., 555.};
      colind = {3, 12};
      A_->insertGlobalValues(it, 2, vals.data(), colind.data());
    }
  }
  A_->fillingIsCompleted();
  //A_->data()->Print(std::cout);

  std::cout << A_->data()->NumGlobalRows() << " " << 
               A_->data()->NumGlobalCols() << std::endl;

  //--------------------------------------------
  //--------------------------------------------
  {
    double o_v[3];
    int o_i[3];
    std::vector<int> nonfilledRows({1,2,4});
    for (auto const & it : nonfilledRows){
      int nn = -1;
      int ec = A_->data()->ExtractGlobalRowCopy( it, 3, nn, o_v, o_i );
      if (ec==0){
	EXPECT_NEAR( 0., o_v[0], 1e-12);
	EXPECT_NEAR( 0., o_v[1], 1e-12);
	EXPECT_NEAR( 0., o_v[2], 1e-12);
      }
      else{
	// this should not work because these rows do not exist
	EXPECT_NE(ec, 0);
      }
    }
  }
  
  // check that row 0 matches
  {
    double o_v[3]; int o_i[3]; int n;
    int ec = A_->data()->ExtractGlobalRowCopy( 0, 3, n, o_v, o_i );
    if (MyPID_ == 0){
      EXPECT_EQ(ec, 0);
      EXPECT_EQ(n, 3);
      EXPECT_DOUBLE_EQ( o_v[0], 1.0);
      EXPECT_DOUBLE_EQ( o_v[1], 5.0);
      EXPECT_DOUBLE_EQ( o_v[2], 2.0);
      EXPECT_EQ( o_i[0], 0);
      EXPECT_EQ( o_i[1], 3);
      EXPECT_EQ( o_i[2], 4);
    }
    else
      EXPECT_EQ(ec, -1);
  }

  // check that row 3 matches
  {
    double o_v[3]; int o_i[3]; int n;
    int ec = A_->data()->ExtractGlobalRowCopy( 3, 3, n, o_v, o_i );
    if (MyPID_ == 0){
      EXPECT_EQ(ec, 0);
      EXPECT_EQ(n, 3);
      EXPECT_DOUBLE_EQ( o_v[0], 11.0);
      EXPECT_DOUBLE_EQ( o_v[1], 55.0);
      EXPECT_DOUBLE_EQ( o_v[2], 22.0);
      EXPECT_EQ( o_i[0], 6);
      EXPECT_EQ( o_i[1], 8);
      EXPECT_EQ( o_i[2], 11);
    }
    else
      EXPECT_EQ(ec, -1);
  }

  // check that row 11 matches
  {
    double o_v[2]; int o_i[2]; int n;
    int ec = A_->data()->ExtractGlobalRowCopy( 11, 2, n, o_v, o_i );
    if (MyPID_ == 2){
      EXPECT_EQ(ec, 0);
      EXPECT_EQ(n , 2);
      EXPECT_DOUBLE_EQ( o_v[0], 555.0);
      EXPECT_DOUBLE_EQ( o_v[1], 111.0);
      EXPECT_EQ( o_i[0], 12);
      EXPECT_EQ( o_i[1], 3);
    }
    else
      EXPECT_EQ(ec, -1);
  }
  
  //--------------------------------------------
  //--------------------------------------------

  // -------------
  // -------------
  // TRANSPOSE
  // TRANSPOSE
  // -------------
  // -------------
  auto At = algebra::mat_ops::transpose(*A_);
  //  At.data()->Print(std::cout);

  // check that row 0 matches
  {
    double o_v[1]; int o_i[1]; int n;
    int ec = At.data()->ExtractGlobalRowCopy( 0, 1, n, o_v, o_i );
    if (MyPID_ == 0){
      EXPECT_EQ(ec, 0); EXPECT_EQ(n, 1);
      EXPECT_DOUBLE_EQ( o_v[0], 1.0);
      EXPECT_EQ( o_i[0], 0);
    }
    else
      EXPECT_EQ(ec, -1);
  }

  // check that row 3 matches
  {
    double o_v[2]; int o_i[2]; int n;
    int ec = At.data()->ExtractGlobalRowCopy( 3, 2, n, o_v, o_i );
    if (MyPID_ == 0){
      EXPECT_EQ(ec, 0); EXPECT_EQ(n, 2);
      EXPECT_DOUBLE_EQ( o_v[0], 5.0);
      EXPECT_DOUBLE_EQ( o_v[1], 111.0);
      EXPECT_EQ( o_i[0], 0);
      EXPECT_EQ( o_i[1], 11);
    }
    else
      EXPECT_EQ(ec, -1);
  }

  // check that row 4 matches
  {
    double o_v[1]; int o_i[1]; int n;
    int ec = At.data()->ExtractGlobalRowCopy( 4, 1, n, o_v, o_i );
    if (MyPID_ == 0){
      EXPECT_EQ(ec, 0); EXPECT_EQ(n, 1);
      EXPECT_DOUBLE_EQ( o_v[0], 2.0);
      EXPECT_EQ( o_i[0], 0);
    }
    else
      EXPECT_EQ(ec, -1);
  }

  // check that row 6 matches
  {
    double o_v[1]; int o_i[1]; int n;
    int ec = At.data()->ExtractGlobalRowCopy( 6, 1, n, o_v, o_i );
    if (MyPID_ == 1){
      EXPECT_EQ(ec, 0); EXPECT_EQ(n, 1);
      EXPECT_DOUBLE_EQ( o_v[0], 11.0);
      EXPECT_EQ( o_i[0], 3);
    }
    else
      EXPECT_EQ(ec, -1);
  }

  // check that row 8 matches
  {
    double o_v[1]; int o_i[1]; int n;
    int ec = At.data()->ExtractGlobalRowCopy( 8, 1, n, o_v, o_i );
    if (MyPID_ == 1){
      EXPECT_EQ(ec, 0); EXPECT_EQ(n, 1);
      EXPECT_DOUBLE_EQ( o_v[0], 55.0);
      EXPECT_EQ( o_i[0], 3);
    }
    else
      EXPECT_EQ(ec, -1);
  }

  // check that row 11 matches
  {
    double o_v[1]; int o_i[1]; int n;
    int ec = At.data()->ExtractGlobalRowCopy( 11, 1, n, o_v, o_i );
    if (MyPID_ == 2){
      EXPECT_EQ(ec, 0); EXPECT_EQ(n, 1);
      EXPECT_DOUBLE_EQ( o_v[0], 22.0);
      EXPECT_EQ( o_i[0], 3);
    }
    else
      EXPECT_EQ(ec, -1);
  }
  
}

