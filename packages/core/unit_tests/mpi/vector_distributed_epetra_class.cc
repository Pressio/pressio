
#include <gtest/gtest.h>
#include "Epetra_MpiComm.h"
#include "CORE_VECTOR"

struct core_vector_distributed_epetraFix
  : public ::testing::Test{
public:
  int rank_;
  Epetra_MpiComm * Comm_;
  int MyPID_;
  int NumProc_;
  const int localSize_ = 5;
  int numGlobalEntries_;
  Epetra_Map * contigMap_;
  Epetra_Vector * x_;
  
  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    Comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
    MyPID_ = Comm_->MyPID();
    NumProc_ = Comm_->NumProc();
    numGlobalEntries_ = Comm_->NumProc() * localSize_;
    contigMap_ = new Epetra_Map(numGlobalEntries_, 0, *Comm_);
    x_ = new Epetra_Vector(*contigMap_);
  }
  
  virtual void TearDown(){
    delete Comm_;
    delete contigMap_;
    delete x_;
  }
};



TEST_F(core_vector_distributed_epetraFix, Constructor){
  using namespace rompp;

  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t a( *contigMap_ );
  ASSERT_EQ( a.globalSize(), numGlobalEntries_ );
  ASSERT_EQ( a.localSize(), localSize_ );
  myvec_t a2( *x_ );
  ASSERT_EQ( a2.globalSize(), numGlobalEntries_ );
  ASSERT_EQ( a2.localSize(), localSize_ );
  myvec_t a3( a2 );
  ASSERT_EQ( a3.globalSize(), numGlobalEntries_ );
  ASSERT_EQ( a3.localSize(), localSize_ );
}


TEST_F(core_vector_distributed_epetraFix,
       QueryWrappedData)
{
  using namespace rompp;

  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *contigMap_ );
  ::testing::StaticAssertTypeEq<decltype(v1.data()),
  				Epetra_Vector * >(); 
  const myvec_t v2( *x_ );
  ::testing::StaticAssertTypeEq< decltype(v2.data()),
  				 const Epetra_Vector * >();
}

TEST_F(core_vector_distributed_epetraFix,
       SubscriptOperator)
{
  using namespace rompp;

  using myvec_t = core::Vector<Epetra_Vector>;

  x_->PutScalar(11.2);
  myvec_t v1( *x_ );
  for (int i=0; i<v1.localSize(); i++){
    v1[i] = 11.2;
  }
  for (int i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( v1[i], (*x_)[i] );
  }
  v1[3] = 56.;
  EXPECT_DOUBLE_EQ( v1[3], 56.0);
}

TEST_F(core_vector_distributed_epetraFix,
       SetScalar){
  using namespace rompp;

  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *contigMap_ );
  v1.putScalar(43.3);

  for (int i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( v1[i], 43.3 );
  }
}

TEST_F(core_vector_distributed_epetraFix,
       AdditionOperator){
  using namespace rompp;

  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *contigMap_ );
  double rankD = static_cast<double>(rank_);
  v1.putScalar( 3.3 +rankD );
  
  myvec_t v2( *contigMap_ );
  v2.putScalar(1.0);

  myvec_t v3 = v1 + v2;
  for (int i=0; i<v3.localSize(); i++){
    EXPECT_DOUBLE_EQ( v3[i], 4.3 + rankD );
  }
  //missing test for a case where vectors are incompatible
}

TEST_F(core_vector_distributed_epetraFix,
       SubtractOperator){
  using namespace rompp;

  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *contigMap_ );
  double rankD = static_cast<double>(rank_);  
  v1.putScalar( 3.3 +rankD );  

  myvec_t v2( *contigMap_ );
  v2.putScalar(1.0);

  myvec_t v3 = v1 - v2;
  for (int i=0; i<v3.localSize(); i++){
    EXPECT_DOUBLE_EQ( v3[i], 2.3 + rankD );
  }
  //missing test for a case where vectors are incompatible
}


TEST_F(core_vector_distributed_epetraFix,
       StarOperator){
  using namespace rompp;

  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *contigMap_ );
  double rankD = static_cast<double>(rank_);  
  v1.putScalar( 3. +rankD );
  
  myvec_t v2( *contigMap_ );
  v2.putScalar(1.0);

  myvec_t v3 = v1 * v2;
  for (int i=0; i<v3.localSize(); i++){
    if (rank_==0)
      EXPECT_DOUBLE_EQ( v3[i], 3. );
    if (rank_==1)
      EXPECT_DOUBLE_EQ( v3[i], 4. );
    if (rank_==2)
      EXPECT_DOUBLE_EQ( v3[i], 5. );
  }
  //missing test for a case where vectors are incompatible
}


TEST_F(core_vector_distributed_epetraFix,
       CompoundAssignAddOperator){
  using namespace rompp;

  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *contigMap_ );
  v1.putScalar( 3. );
  
  myvec_t v2( *contigMap_ );
  v2.putScalar(1.0);

  v1 += v2;
  for (int i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( v1[i], 4. );
  }
  //missing test for a case where vectors are incompatible
}


TEST_F(core_vector_distributed_epetraFix,
       CompoundAssignSubtractOperator){
  using namespace rompp;

  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *contigMap_ );
  v1.putScalar( 3. );
  
  myvec_t v2( *contigMap_ );
  v2.putScalar(1.0);

  v1 -= v2;
  for (int i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( v1[i], 2. );
  }
  //missing test for a case where vectors are incompatible
}


TEST_F(core_vector_distributed_epetraFix,
       SetZero){
  using namespace rompp;

  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *contigMap_ );
  v1.setZero();
  
  for (int i=0; i<v1.localSize(); i++){
    EXPECT_NEAR( v1[i], 0.0, 1e-12 );
  }
}


TEST_F(core_vector_distributed_epetraFix,
       Empty){
  using namespace rompp;

  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *contigMap_ );
  ASSERT_FALSE( v1.empty() );
}


TEST_F(core_vector_distributed_epetraFix,
       replaceGlobalData){

  using namespace rompp;

  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *contigMap_ );
  v1.setZero();

  int id[2];
  double val[2];
  if(rank_==0){
    id[0] = 0; val[0]=1.2;
    id[1] = 2; val[1]=3.1;
    v1.replaceGlobalValues(2, id, val);

    EXPECT_DOUBLE_EQ( v1[0], 1.2 );
    EXPECT_DOUBLE_EQ( v1[1], 0.0 );
    EXPECT_DOUBLE_EQ( v1[2], 3.1 );
  }
}


TEST_F(core_vector_distributed_epetraFix,
       InPlaceOp){
  using namespace rompp;

  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *contigMap_ );
  v1.putScalar(1.1);

  myvec_t v2( *contigMap_ );
  v2.putScalar(2.1);
  
  v1.inPlaceOp<std::plus<double>>(2.0, 0.0, v2);

  EXPECT_DOUBLE_EQ( v1[0], 2.2);
  EXPECT_DOUBLE_EQ( v1[1], 2.2);
  EXPECT_DOUBLE_EQ( v1[2], 2.2);
  EXPECT_DOUBLE_EQ( v1[3], 2.2);
  EXPECT_DOUBLE_EQ( v1[4], 2.2);
}


TEST_F(core_vector_distributed_epetraFix,
       Scale){
  using namespace rompp;

  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *contigMap_ );
  v1.putScalar(1.1);
  v1.scale(3.);

  EXPECT_DOUBLE_EQ( v1[0], 3.3);
  EXPECT_DOUBLE_EQ( v1[1], 3.3);
  EXPECT_DOUBLE_EQ( v1[2], 3.3);
  EXPECT_DOUBLE_EQ( v1[3], 3.3);
  EXPECT_DOUBLE_EQ( v1[4], 3.3);
}









// using core_vector_distributed_epetra_DeathTest
// = core_vector_distributed_epetraFix;
// TEST_F(core_vector_distributed_epetra_DeathTest,
//        EpetraVectorSubscriptOperator)
// {
//   if (rank_==0){
//     using myvec_t = core::Vector<Epetra_Vector>;
//     myvec_t v1( *getMap() );
//     int localSize = numLocalEntries();
//     ASSERT_DEATH(v1[localSize+1]==4.0, "Exit code:    1");
//   }
// }
