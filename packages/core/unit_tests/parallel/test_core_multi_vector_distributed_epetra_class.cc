
#include <gtest/gtest.h>
#include "Epetra_MpiComm.h"
#include "CORE_VECTOR"
#include "CORE_MULTI_VECTOR"

struct core_multi_vector_distributed_epetraFix
  : public ::testing::Test{
public:
  int rank_;
  Epetra_MpiComm * Comm_;
  int MyPID_;
  int NumProc_;
  const int localSize_ = 5;
  const int numVectors_ = 3;
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

TEST_F(core_multi_vector_distributed_epetraFix,
       EpetraMultiVectorConstructor)
{
  using nat_t = Epetra_MultiVector;
  using mymvec_t = core::MultiVector<nat_t>;
  STATIC_ASSERT_IS_CORE_MULTI_VECTOR_WRAPPER(mymvec_t);

  mymvec_t a( *getMap(), numVectors_ );
  ASSERT_FALSE( a.empty() );

  mymvec_t b( getNative() );  
  ASSERT_FALSE( b.empty() );

  EXPECT_EQ( b.globalNumVectors(), 3 );
  EXPECT_EQ( b.localNumVectors(), 3 );
  EXPECT_EQ( b.globalLength(), 15 );
  EXPECT_EQ( b.localLength(), 5);

  for (int i=0; i<localSize_; i++)
    for (int j=0; j<b.globalNumVectors(); j++)
      EXPECT_NEAR( 0.0, b(i,j), 1e-12);

  if(MyPID_==0)
    b.replaceGlobalValue(1,1, 43.3);
  if(MyPID_==0)
    EXPECT_NEAR( 43.3, b(1,1), 1e-12);
  else
    EXPECT_NEAR( 0.0, b(1,1), 1e-12);

  b.scale(2.0);
  if(MyPID_==0)
    EXPECT_NEAR( 86.6, b(1,1), 1e-12);
  
}




// TEST_F(core_vector_distributed_epetraFix,
//        EpetraVectorQueryWrappedData)
// {    
//   using myvec_t = core::Vector<Epetra_Vector>;
//   myvec_t v1( *getMap() );
//   ::testing::StaticAssertTypeEq<decltype(v1.data()),
//   				Epetra_Vector * >(); 
//   const myvec_t v2( *getVector() );
//   ::testing::StaticAssertTypeEq< decltype(v2.data()),
//   				 const Epetra_Vector * >();
// }

// TEST_F(core_vector_distributed_epetraFix,
//        EpetraVectorSubscriptOperator)
// {
//   using myvec_t = core::Vector<Epetra_Vector>;

//   fillScalar(11.2);
//   myvec_t v1( *getMap() );
//   for (int i=0; i<v1.localSize(); i++){
//     v1[i] = 11.2;
//   }
//   for (int i=0; i<v1.localSize(); i++){
//     EXPECT_DOUBLE_EQ( v1[i], (*getVector())[i] );
//   }
//   v1[3] = 56.;
//   EXPECT_DOUBLE_EQ( v1[3], 56.0);
// }


// TEST_F(core_vector_distributed_epetraFix,
//        EpetraVectorSetScalar)
// {
//   using myvec_t = core::Vector<Epetra_Vector>;
//   myvec_t v1( *getMap() );
//   v1.putScalar(43.3);

//   for (int i=0; i<v1.localSize(); i++){
//     EXPECT_DOUBLE_EQ( v1[i], 43.3 );
//   }
// }


// TEST_F(core_vector_distributed_epetraFix,
//        EpetraVectorAdditionOperator)
// {
//   using myvec_t = core::Vector<Epetra_Vector>;
//   myvec_t v1( *getMap() );
//   double rankD = static_cast<double>(getRank());  
//   v1.putScalar( 3.3 +rankD );
  
//   myvec_t v2( *getMap() );
//   v2.putScalar(1.0);

//   myvec_t v3 = v1 + v2;
//   for (int i=0; i<v3.localSize(); i++){
//     EXPECT_DOUBLE_EQ( v3[i], 4.3 + rankD );
//   }
//   //missing test for a case where vectors are incompatible
// }


// TEST_F(core_vector_distributed_epetraFix,
//        EpetraVectorSubtractOperator)
// {
//   using myvec_t = core::Vector<Epetra_Vector>;
//   myvec_t v1( *getMap() );
//   double rankD = static_cast<double>(getRank());  
//   v1.putScalar( 3.3 +rankD );
  
//   myvec_t v2( *getMap() );
//   v2.putScalar(1.0);

//   myvec_t v3 = v1 - v2;
//   for (int i=0; i<v3.localSize(); i++){
//     EXPECT_DOUBLE_EQ( v3[i], 2.3 + rankD );
//   }
//   //missing test for a case where vectors are incompatible
// }


// TEST_F(core_vector_distributed_epetraFix,
//        EpetraVectorStarOperator)
// {
//   using myvec_t = core::Vector<Epetra_Vector>;
//   myvec_t v1( *getMap() );
//   double rankD = static_cast<double>(getRank());  
//   v1.putScalar( 3. +rankD );
  
//   myvec_t v2( *getMap() );
//   v2.putScalar(1.0);

//   myvec_t v3 = v1 * v2;
//   for (int i=0; i<v3.localSize(); i++){
//     if (getRank()==0)
//       EXPECT_DOUBLE_EQ( v3[i], 3. );
//     if (getRank()==1)
//       EXPECT_DOUBLE_EQ( v3[i], 4. );
//     if (getRank()==2)
//       EXPECT_DOUBLE_EQ( v3[i], 5. );
//   }
//   //missing test for a case where vectors are incompatible
// }


// TEST_F(core_vector_distributed_epetraFix,
//        EpetraVectorCompoundAssignAddOperator)
// {
//   using myvec_t = core::Vector<Epetra_Vector>;
//   myvec_t v1( *getMap() );
//   v1.putScalar( 3. );
  
//   myvec_t v2( *getMap() );
//   v2.putScalar(1.0);

//   v1 += v2;
//   for (int i=0; i<v1.localSize(); i++){
//     EXPECT_DOUBLE_EQ( v1[i], 4. );
//   }
//   //missing test for a case where vectors are incompatible
// }


// TEST_F(core_vector_distributed_epetraFix,
//        EpetraVectorCompoundAssignSubtractOperator)
// {
//   using myvec_t = core::Vector<Epetra_Vector>;
//   myvec_t v1( *getMap() );
//   v1.putScalar( 3. );
  
//   myvec_t v2( *getMap() );
//   v2.putScalar(1.0);

//   v1 -= v2;
//   for (int i=0; i<v1.localSize(); i++){
//     EXPECT_DOUBLE_EQ( v1[i], 2. );
//   }
//   //missing test for a case where vectors are incompatible
// }

