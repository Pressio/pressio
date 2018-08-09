
#include <gtest/gtest.h>
#include "Epetra_MpiComm.h"
#include "CORE_VECTOR"
#include "CORE_MATRIX"

struct core_matrix_sparse_distributed_epetraFix
  : public ::testing::Test{
public:
  int rank_;
  Epetra_MpiComm * Comm_;
  int MyPID_;
  int NumProc_;
  const int localSize_ = 2;
  int numGlobalEntries_;
  Epetra_Map * contigMap_;
  core::Matrix<Epetra_CrsMatrix> * A_;
  core::Matrix<Epetra_CrsMatrix> * B_;
  
  virtual void SetUp()
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    Comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
    MyPID_ = Comm_->MyPID();
    NumProc_ = Comm_->NumProc();
    numGlobalEntries_ = Comm_->NumProc() * localSize_;
    contigMap_ = new Epetra_Map(numGlobalEntries_, 0, *Comm_);
    A_ = new core::Matrix<Epetra_CrsMatrix>(*contigMap_, 3);
    B_ = new core::Matrix<Epetra_CrsMatrix>(*contigMap_, 4);
  }
  
  virtual void TearDown(){
    delete Comm_;
    delete contigMap_;
    delete A_;
    delete B_;
  }
};

TEST_F(core_matrix_sparse_distributed_epetraFix,Test1)
{
  //-----------
  // FILL A 
  //-----------
  {
    int myN = contigMap_->NumMyElements();
    // get my global indices
    std::vector<int> mygid(myN);
    contigMap_->MyGlobalElements( mygid.data() );

    std::array<double,3> vals;
    std::array<int,3> colind;
    for (auto const & it : mygid){
      if(it == 0){
	vals = {1., 2.};
	colind = {0, 2};
	A_->insertGlobalValues(it, 2, vals.data(), colind.data());
      }
      if(it == 2){
	vals = {1.};
	colind = {2};
	A_->insertGlobalValues(it, 1, vals.data(), colind.data());
      }
    }
    A_->fillingIsCompleted();
    //    A_->data()->Print(std::cout);
  }

  //-----------
  // FILL B
  //-----------
  {
    int myN = contigMap_->NumMyElements();
    // get my global indices
    std::vector<int> mygid(myN);
    contigMap_->MyGlobalElements( mygid.data() );
    std::array<double,4> vals;
    std::array<int,4> colind;
    for (auto const & it : mygid){
      if(it == 0){
	vals = {2.}; colind = {0};
	B_->insertGlobalValues(it, 1, vals.data(), colind.data());
      }
      if(it == 2){
	vals = {1., 6.}; colind = {0, 2};
	B_->insertGlobalValues(it, 2, vals.data(), colind.data());
      }
      if(it == 3){
	vals = {4.}; colind = {0};
	B_->insertGlobalValues(it, 1, vals.data(), colind.data());
      }
    }
    B_->fillingIsCompleted();
    //    B_->data()->Print(std::cout);
    //    std::cout << "GIGI " << B_->globalRows() << std::endl;
  }

  //-----------
  // product 1
  //-----------
  {
    core::Matrix<Epetra_CrsMatrix> myC_(*contigMap_, 3);
    matrixMatrixProduct(*A_, *B_, myC_, false, false);
    //myC_.data()->Print(std::cout);
    {
      double o_v[2]; int o_i[2]; int n;
      int ec = myC_.data()->ExtractGlobalRowCopy( 0, 2, n, o_v, o_i );
      if (MyPID_ == 0){
	EXPECT_EQ(ec, 0); EXPECT_EQ(n, 2);
	EXPECT_DOUBLE_EQ( o_v[0], 4.0);
	EXPECT_DOUBLE_EQ( o_v[1], 12.0);
	EXPECT_EQ( o_i[0], 0);
	EXPECT_EQ( o_i[1], 2);
      }
      else
	EXPECT_EQ(ec, -1);
    }

    {
      double o_v[2]; int o_i[2]; int n;
      int ec = myC_.data()->ExtractGlobalRowCopy( 2, 2, n, o_v, o_i );
      if (MyPID_ == 1){
	EXPECT_EQ(ec, 0); EXPECT_EQ(n, 2);
	EXPECT_DOUBLE_EQ( o_v[0], 6.0);
	EXPECT_DOUBLE_EQ( o_v[1], 1.0);
	EXPECT_EQ( o_i[0], 2);
	EXPECT_EQ( o_i[1], 0);
      }
      else
	EXPECT_EQ(ec, -1);
    }
  }

  //-----------
  // product 2
  //-----------
  {
    auto myC_ = matrixMatrixProduct(*A_, *B_, false, false);
    {
      double o_v[2]; int o_i[2]; int n;
      int ec = myC_.data()->ExtractGlobalRowCopy( 0, 2, n, o_v, o_i );
      if (MyPID_ == 0){
	EXPECT_EQ(ec, 0); EXPECT_EQ(n, 2);
	EXPECT_DOUBLE_EQ( o_v[0], 4.0);
	EXPECT_DOUBLE_EQ( o_v[1], 12.0);
	EXPECT_EQ( o_i[0], 0);
	EXPECT_EQ( o_i[1], 2);
      }
      else
	EXPECT_EQ(ec, -1);
    }

    {
      double o_v[2]; int o_i[2]; int n;
      int ec = myC_.data()->ExtractGlobalRowCopy( 2, 2, n, o_v, o_i );
      if (MyPID_ == 1){
	EXPECT_EQ(ec, 0); EXPECT_EQ(n, 2);
	EXPECT_DOUBLE_EQ( o_v[0], 6.);
	EXPECT_DOUBLE_EQ( o_v[1], 1.0);
	EXPECT_EQ( o_i[0], 2);
	EXPECT_EQ( o_i[1], 0);
      }
      else
	EXPECT_EQ(ec, -1);
    }
  }
    
}
