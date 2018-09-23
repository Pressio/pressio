
#include <gtest/gtest.h>
#include "Epetra_MpiComm.h"
#include "CORE_ALL"

struct epetraFix
  : public ::testing::Test{
public:
  int rank_;
  Epetra_MpiComm * Comm_;
  int MyPID_;
  int NumProc_;
  const int localSize_ = 3;
  const int numVectors_ = 4;
  int numGlobalEntries_;
  Epetra_Map * dataMap_;
  Epetra_MultiVector * mv_;
  Epetra_Vector * x_;
  
  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    Comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
    MyPID_ = Comm_->MyPID();
    NumProc_ = Comm_->NumProc();
    numGlobalEntries_ = Comm_->NumProc() * localSize_;
    dataMap_ = new Epetra_Map(numGlobalEntries_, 0, *Comm_);
    mv_ = new Epetra_MultiVector(*dataMap_, numVectors_);
    x_ = new Epetra_Vector(*dataMap_);
  }

  virtual void TearDown(){
    delete Comm_;
    delete dataMap_;
    delete mv_;
    delete x_;
  }
};

TEST_F(epetraFix, MVVecDotProduct){

  using namespace rompp;

  assert(NumProc_ == 3);
  using mvec_t = core::MultiVector<Epetra_MultiVector>;
  STATIC_ASSERT_IS_CORE_MULTI_VECTOR_WRAPPER(mvec_t);
  using vec_t = core::Vector<Epetra_Vector>;
  STATIC_ASSERT_IS_CORE_VECTOR_WRAPPER(vec_t);
  mvec_t MV(*mv_);
  vec_t b(*x_);

  EXPECT_EQ( MV.globalNumVectors(), 4 );
  EXPECT_EQ( MV.localNumVectors(), 4 );
  EXPECT_EQ( MV.globalLength(), 9 );
  EXPECT_EQ( MV.localLength(), 3);

  for (int i=0; i<localSize_; i++)
    for (int j=0; j<MV.globalNumVectors(); j++)
      EXPECT_NEAR( 0.0, MV(i,j), 1e-12);

  if(MyPID_==0){
    MV(0,0) = 3.2;
    MV(1,0) = 1.2;
    MV(2,1) = 4;
    MV(0,1) = 1.2;
  }

  if(MyPID_==1){
    MV(2,2) = 3;
  }
  
  b[0] = 1.0;
  b[1] = 1.0;
  b[2] = 1.0;

  // MV.data()->Print(std::cout);
  // b.data()->Print(std::cout);
  auto res = core::ops::dot(MV, b);

  EXPECT_EQ((int) res.size(), 4);
  EXPECT_NEAR(res[0], 4.4, 1e-12);
  EXPECT_NEAR(res[1], 5.2, 1e-12);
  EXPECT_NEAR(res[2], 3., 1e-12);
  EXPECT_NEAR(res[3], 0., 1e-12);

  // if (MyPID_==0){
  //   //std::cout << res << " ";
  //   for (const auto & it : res)
  //     std::cout << it << " ";
  //   std::cout << "\n";
  // }  
}
