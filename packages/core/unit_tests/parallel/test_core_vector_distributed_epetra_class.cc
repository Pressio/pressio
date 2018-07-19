
#include <gtest/gtest.h>
// #include "vector/meta/core_vector_meta.hpp"
// #include "vector/concrete/core_vector_distributed_epetra.hpp"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"
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
  
  virtual void SetUp()
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    Comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
    MyPID_ = Comm_->MyPID();
    NumProc_ = Comm_->NumProc();
    numGlobalEntries_ = Comm_->NumProc() * localSize_;
    contigMap_ = new Epetra_Map(numGlobalEntries_, 0, *Comm_);
    x_ = new Epetra_Vector(*contigMap_);
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
    return contigMap_;
  }    
  const Epetra_Vector * getVector() const{
    return x_;
  }
  void fillScalar(double value){
    x_->PutScalar(value);
  }
    
  virtual void TearDown(){
    delete Comm_;
    delete contigMap_;
    delete x_;
  }
};



using core_vector_distributed_epetra_DeathTest
= core_vector_distributed_epetraFix;

TEST_F(core_vector_distributed_epetra_DeathTest,
       EpetraVectorSubscriptOperator)
{
  if (getRank()==0){
    using myvec_t = core::Vector<Epetra_Vector>;
    myvec_t v1( *getMap() );
    int localSize = numLocalEntries();
    ASSERT_DEATH(v1[localSize+1]=4.0, "Assertion failed:");
  }
}



TEST_F(core_vector_distributed_epetraFix,
       EpetraVectorConstructor)
{
  using myvec_t = core::Vector<Epetra_Vector>;
  ASSERT_TRUE( core::details::traits<myvec_t>::isEpetra==1 );

  myvec_t a( *getMap() );
  ASSERT_EQ( a.globalSize(), numGlobalEntries() );
  ASSERT_EQ( a.localSize(), numLocalEntries() );
  myvec_t a2( *getVector() );
  ASSERT_EQ( a2.globalSize(), numGlobalEntries() );
  ASSERT_EQ( a2.localSize(), numLocalEntries() );
}


TEST_F(core_vector_distributed_epetraFix,
       EpetraVectorQueryWrappedData)
{    
  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *getMap() );
  ::testing::StaticAssertTypeEq<decltype(v1.data()),
  				Epetra_Vector * >(); 
  const myvec_t v2( *getVector() );
  ::testing::StaticAssertTypeEq< decltype(v2.data()),
  				 const Epetra_Vector * >();
}

TEST_F(core_vector_distributed_epetraFix,
       EpetraVectorSubscriptOperator)
{
  using myvec_t = core::Vector<Epetra_Vector>;

  fillScalar(11.2);
  myvec_t v1( *getMap() );
  for (int i=0; i<v1.localSize(); i++){
    v1[i] = 11.2;
  }
  for (int i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( v1[i], (*getVector())[i] );
  }
  v1[3] = 56.;
  EXPECT_DOUBLE_EQ( v1[3], 56.0);
}


TEST_F(core_vector_distributed_epetraFix,
       EpetraVectorSetScalar)
{
  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *getMap() );
  v1.putScalar(43.3);

  for (int i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( v1[i], 43.3 );
  }
}


TEST_F(core_vector_distributed_epetraFix,
       EpetraVectorAdditionOperator)
{
  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *getMap() );
  double rankD = static_cast<double>(getRank());  
  v1.putScalar( 3.3 +rankD );
  
  myvec_t v2( *getMap() );
  v2.putScalar(1.0);

  myvec_t v3 = v1 + v2;
  for (int i=0; i<v3.localSize(); i++){
    EXPECT_DOUBLE_EQ( v3[i], 4.3 + rankD );
  }
  //missing test for a case where vectors are incompatible
}


TEST_F(core_vector_distributed_epetraFix,
       EpetraVectorSubtractOperator)
{
  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *getMap() );
  double rankD = static_cast<double>(getRank());  
  v1.putScalar( 3.3 +rankD );
  
  myvec_t v2( *getMap() );
  v2.putScalar(1.0);

  myvec_t v3 = v1 - v2;
  for (int i=0; i<v3.localSize(); i++){
    EXPECT_DOUBLE_EQ( v3[i], 2.3 + rankD );
  }
  //missing test for a case where vectors are incompatible
}


TEST_F(core_vector_distributed_epetraFix,
       EpetraVectorStarOperator)
{
  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *getMap() );
  double rankD = static_cast<double>(getRank());  
  v1.putScalar( 3. +rankD );
  
  myvec_t v2( *getMap() );
  v2.putScalar(1.0);

  myvec_t v3 = v1 * v2;
  for (int i=0; i<v3.localSize(); i++){
    if (getRank()==0)
      EXPECT_DOUBLE_EQ( v3[i], 3. );
    if (getRank()==1)
      EXPECT_DOUBLE_EQ( v3[i], 4. );
    if (getRank()==2)
      EXPECT_DOUBLE_EQ( v3[i], 5. );
  }
  //missing test for a case where vectors are incompatible
}


TEST_F(core_vector_distributed_epetraFix,
       EpetraVectorCompoundAssignAddOperator)
{
  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *getMap() );
  v1.putScalar( 3. );
  
  myvec_t v2( *getMap() );
  v2.putScalar(1.0);

  v1 += v2;
  for (int i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( v1[i], 4. );
  }
  //missing test for a case where vectors are incompatible
}


TEST_F(core_vector_distributed_epetraFix,
       EpetraVectorCompoundAssignSubtractOperator)
{
  using myvec_t = core::Vector<Epetra_Vector>;
  myvec_t v1( *getMap() );
  v1.putScalar( 3. );
  
  myvec_t v2( *getMap() );
  v2.putScalar(1.0);

  v1 -= v2;
  for (int i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( v1[i], 2. );
  }
  //missing test for a case where vectors are incompatible
}

