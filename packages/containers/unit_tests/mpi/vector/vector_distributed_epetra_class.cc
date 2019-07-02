
#include "epetra_only_fixtures.hpp"

TEST_F(epetraVectorGlobSize15Fixture, Constructor){
  using namespace pressio;

  using myvec_t = containers::Vector<Epetra_Vector>;
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


// TEST_F(epetraVectorGlobSize15Fixture,
//        QueryWrappedData)
// {
//   using namespace pressio;

//   using myvec_t = containers::Vector<Epetra_Vector>;
//   myvec_t v1( *contigMap_ );
//   ::testing::StaticAssertTypeEq<decltype(v1.data()),
//   				Epetra_Vector * >(); 
//   const myvec_t v2( *x_ );
//   ::testing::StaticAssertTypeEq< decltype(v2.data()),
//   				 const Epetra_Vector * >();
// }


TEST_F(epetraVectorGlobSize15Fixture,
       SubscriptOperator)
{
  using namespace pressio;

  using myvec_t = containers::Vector<Epetra_Vector>;

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


// TEST_F(epetraVectorGlobSize15Fixture,
//        SetScalar){
//   using namespace pressio;

//   using myvec_t = containers::Vector<Epetra_Vector>;
//   myvec_t v1( *contigMap_ );
//   v1.putScalar(43.3);

//   for (int i=0; i<v1.localSize(); i++){
//     EXPECT_DOUBLE_EQ( v1[i], 43.3 );
//   }
// }

// TEST_F(epetraVectorGlobSize15Fixture,
//        AdditionOperator){
//   using namespace pressio;

//   using myvec_t = containers::Vector<Epetra_Vector>;
//   myvec_t v1( *contigMap_ );
//   double rankD = static_cast<double>(rank_);
//   v1.putScalar( 3.3 +rankD );
  
//   myvec_t v2( *contigMap_ );
//   v2.putScalar(1.0);
  
//   myvec_t v3 = v1 + v2;
//   for (int i=0; i<v3.localSize(); i++){
//     EXPECT_DOUBLE_EQ( v3[i], 4.3 + rankD );
//   }
//   //missing test for a case where vectors are incompatible
// }


TEST_F(epetraVectorGlobSize15Fixture,
       expreTempPlus){
  using namespace pressio;

  using myvec_t = containers::Vector<Epetra_Vector>;
  myvec_t v1( *contigMap_ );
  double rankD = static_cast<double>(rank_);
  v1.putScalar( 3.3 +rankD );
  
  myvec_t v2( *contigMap_ );
  v2.putScalar(1.0);

  myvec_t v3( *contigMap_ );
  v3.putScalar(1.2);

  myvec_t v4( *contigMap_ );

  v4 = (v1+v2)*2.;
  v4.data()->Print(std::cout);

  v4 = 3.*(v1+v2);
  v4.data()->Print(std::cout);
  
  v4 = 3.*v2;
  v4.data()->Print(std::cout);

  v4 = v2*2.;
  v4.data()->Print(std::cout);
  
  v4 = v2 + v1;
  v4.data()->Print(std::cout);

  v4 = v2 - v1;
  v4.data()->Print(std::cout);
  
  v4 = v2 + v1 + v3;
  v4.data()->Print(std::cout);
  
  v4 = v2 + v1 + v3;
  v4.data()->Print(std::cout);

  v4 = v2*2. + 1.*v1;
  v4.data()->Print(std::cout);

  v4 = v2*2. + 1.*v1 + v3;
  v4.data()->Print(std::cout);
  
  // for (int i=0; i<v3.localSize(); i++){
  //   EXPECT_DOUBLE_EQ( v3[i], 4.3 + rankD );
  // }
}




// TEST_F(epetraVectorGlobSize15Fixture,
//        SubtractOperator){
//   using namespace pressio;

//   using myvec_t = containers::Vector<Epetra_Vector>;
//   myvec_t v1( *contigMap_ );
//   double rankD = static_cast<double>(rank_);  
//   v1.putScalar( 3.3 +rankD );  

//   myvec_t v2( *contigMap_ );
//   v2.putScalar(1.0);

//   myvec_t v3 = v1 - v2;
//   for (int i=0; i<v3.localSize(); i++){
//     EXPECT_DOUBLE_EQ( v3[i], 2.3 + rankD );
//   }
//   //missing test for a case where vectors are incompatible
// }


// TEST_F(epetraVectorGlobSize15Fixture,
//        StarOperator){
//   using namespace pressio;

//   using myvec_t = containers::Vector<Epetra_Vector>;
//   myvec_t v1( *contigMap_ );
//   double rankD = static_cast<double>(rank_);  
//   v1.putScalar( 3. +rankD );
  
//   myvec_t v2( *contigMap_ );
//   v2.putScalar(1.0);

//   myvec_t v3 = v1 * v2;
//   for (int i=0; i<v3.localSize(); i++){
//     if (rank_==0)
//       EXPECT_DOUBLE_EQ( v3[i], 3. );
//     if (rank_==1)
//       EXPECT_DOUBLE_EQ( v3[i], 4. );
//     if (rank_==2)
//       EXPECT_DOUBLE_EQ( v3[i], 5. );
//   }
//   //missing test for a case where vectors are incompatible
// }


// TEST_F(epetraVectorGlobSize15Fixture,
//        CompoundAssignAddOperator){
//   using namespace pressio;

//   using myvec_t = containers::Vector<Epetra_Vector>;
//   myvec_t v1( *contigMap_ );
//   v1.putScalar( 3. );
  
//   myvec_t v2( *contigMap_ );
//   v2.putScalar(1.0);

//   v1 += v2;
//   for (int i=0; i<v1.localSize(); i++){
//     EXPECT_DOUBLE_EQ( v1[i], 4. );
//   }
//   //missing test for a case where vectors are incompatible
// }


// TEST_F(epetraVectorGlobSize15Fixture,
//        CompoundAssignSubtractOperator){
//   using namespace pressio;

//   using myvec_t = containers::Vector<Epetra_Vector>;
//   myvec_t v1( *contigMap_ );
//   v1.putScalar( 3. );
  
//   myvec_t v2( *contigMap_ );
//   v2.putScalar(1.0);

//   v1 -= v2;
//   for (int i=0; i<v1.localSize(); i++){
//     EXPECT_DOUBLE_EQ( v1[i], 2. );
//   }
//   //missing test for a case where vectors are incompatible
// }


// TEST_F(epetraVectorGlobSize15Fixture,
//        SetZero){
//   using namespace pressio;

//   using myvec_t = containers::Vector<Epetra_Vector>;
//   myvec_t v1( *contigMap_ );
//   v1.setZero();
  
//   for (int i=0; i<v1.localSize(); i++){
//     EXPECT_NEAR( v1[i], 0.0, 1e-12 );
//   }
// }


// TEST_F(epetraVectorGlobSize15Fixture,
//        Empty){
//   using namespace pressio;

//   using myvec_t = containers::Vector<Epetra_Vector>;
//   myvec_t v1( *contigMap_ );
//   ASSERT_FALSE( v1.empty() );
// }


// TEST_F(epetraVectorGlobSize15Fixture,
//        replaceGlobalData){

//   using namespace pressio;

//   using myvec_t = containers::Vector<Epetra_Vector>;
//   myvec_t v1( *contigMap_ );
//   v1.setZero();

//   int id[2];
//   double val[2];
//   if(rank_==0){
//     id[0] = 0; val[0]=1.2;
//     id[1] = 2; val[1]=3.1;
//     v1.replaceGlobalValues(2, id, val);

//     EXPECT_DOUBLE_EQ( v1[0], 1.2 );
//     EXPECT_DOUBLE_EQ( v1[1], 0.0 );
//     EXPECT_DOUBLE_EQ( v1[2], 3.1 );
//   }
// }


// TEST_F(epetraVectorGlobSize15Fixture,
//        InPlaceOp){
//   using namespace pressio;

//   using myvec_t = containers::Vector<Epetra_Vector>;
//   myvec_t v1( *contigMap_ );
//   v1.putScalar(1.1);

//   myvec_t v2( *contigMap_ );
//   v2.putScalar(2.1);
  
//   v1.inPlaceOp<std::plus<double>>(2.0, 0.0, v2);

//   EXPECT_DOUBLE_EQ( v1[0], 2.2);
//   EXPECT_DOUBLE_EQ( v1[1], 2.2);
//   EXPECT_DOUBLE_EQ( v1[2], 2.2);
//   EXPECT_DOUBLE_EQ( v1[3], 2.2);
//   EXPECT_DOUBLE_EQ( v1[4], 2.2);
// }


// TEST_F(epetraVectorGlobSize15Fixture,
//        Scale){
//   using namespace pressio;

//   using myvec_t = containers::Vector<Epetra_Vector>;
//   myvec_t v1( *contigMap_ );
//   v1.putScalar(1.1);
//   v1.scale(3.);

//   EXPECT_DOUBLE_EQ( v1[0], 3.3);
//   EXPECT_DOUBLE_EQ( v1[1], 3.3);
//   EXPECT_DOUBLE_EQ( v1[2], 3.3);
//   EXPECT_DOUBLE_EQ( v1[3], 3.3);
//   EXPECT_DOUBLE_EQ( v1[4], 3.3);
// }









// using containers_vector_distributed_epetra_DeathTest
// = epetraVectorGlobSize15Fixture;
// TEST_F(containers_vector_distributed_epetra_DeathTest,
//        EpetraVectorSubscriptOperator)
// {
//   if (rank_==0){
//     using myvec_t = containers::Vector<Epetra_Vector>;
//     myvec_t v1( *getMap() );
//     int localSize = numLocalEntries();
//     ASSERT_DEATH(v1[localSize+1]==4.0, "Exit code:    1");
//   }
// }
