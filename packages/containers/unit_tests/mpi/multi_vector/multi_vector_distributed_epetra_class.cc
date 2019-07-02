
#include "epetra_only_fixtures.hpp"


TEST_F(epetraMultiVectorR9C4VecS9Fixture,
       EpetraMultiVectorConstructor){
  
  using namespace pressio;

  using nat_t = Epetra_MultiVector;
  using mymvec_t = containers::MultiVector<nat_t>;
  STATIC_ASSERT_IS_CONTAINERS_MULTI_VECTOR_WRAPPER(mymvec_t);

  mymvec_t a( *dataMap_, numVectors_ );
  ASSERT_FALSE( a.empty() );
  //a.data()->Print(std::cout);

  mymvec_t b( *mv_ );
  ASSERT_FALSE( b.empty() );
  //b.data()->Print(std::cout);

  EXPECT_EQ( b.globalNumVectors(), 4 );
  EXPECT_EQ( b.localNumVectors(), 4 );
  EXPECT_EQ( b.globalLength(), 9 );
  EXPECT_EQ( b.localLength(), 3);

  for (int i=0; i<localSize_; i++)
    for (int j=0; j<b.globalNumVectors(); j++)
      EXPECT_NEAR( 0.0, b(i,j), 1e-12);

  if(rank_==0)
    b.replaceGlobalValue(1,1, 43.3);
  if(rank_==1)
    b.replaceGlobalValue(4,2, 13.3);
  b.data()->Print(std::cout);

  if(rank_==0)
    EXPECT_NEAR( 43.3, b(1,1), 1e-12);
  else if (rank_==1)
    EXPECT_NEAR( 13.3, b(1,2), 1e-12);
  else
    EXPECT_NEAR( 0.0, b(1,1), 1e-12);

  b.scale(2.0);
  if(rank_==0){
    EXPECT_NEAR( 86.6, b(1,1), 1e-12);
  }
  if (rank_==1){
    EXPECT_NEAR( 26.6, b(1,2), 1e-12);
  }
}




// TEST_F(containers_vector_distributed_epetraFix,
//        EpetraVectorQueryWrappedData)
// {    
//   using myvec_t = containers::Vector<Epetra_Vector>;
//   myvec_t v1( *getMap() );
//   ::testing::StaticAssertTypeEq<decltype(v1.data()),
//   				Epetra_Vector * >(); 
//   const myvec_t v2( *getVector() );
//   ::testing::StaticAssertTypeEq< decltype(v2.data()),
//   				 const Epetra_Vector * >();
// }

// TEST_F(containers_vector_distributed_epetraFix,
//        EpetraVectorSubscriptOperator)
// {
//   using myvec_t = containers::Vector<Epetra_Vector>;

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


// TEST_F(containers_vector_distributed_epetraFix,
//        EpetraVectorSetScalar)
// {
//   using myvec_t = containers::Vector<Epetra_Vector>;
//   myvec_t v1( *getMap() );
//   v1.putScalar(43.3);

//   for (int i=0; i<v1.localSize(); i++){
//     EXPECT_DOUBLE_EQ( v1[i], 43.3 );
//   }
// }


// TEST_F(containers_vector_distributed_epetraFix,
//        EpetraVectorAdditionOperator)
// {
//   using myvec_t = containers::Vector<Epetra_Vector>;
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


// TEST_F(containers_vector_distributed_epetraFix,
//        EpetraVectorSubtractOperator)
// {
//   using myvec_t = containers::Vector<Epetra_Vector>;
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


// TEST_F(containers_vector_distributed_epetraFix,
//        EpetraVectorStarOperator)
// {
//   using myvec_t = containers::Vector<Epetra_Vector>;
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


// TEST_F(containers_vector_distributed_epetraFix,
//        EpetraVectorCompoundAssignAddOperator)
// {
//   using myvec_t = containers::Vector<Epetra_Vector>;
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


// TEST_F(containers_vector_distributed_epetraFix,
//        EpetraVectorCompoundAssignSubtractOperator)
// {
//   using myvec_t = containers::Vector<Epetra_Vector>;
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

