
#include "epetra_only_fixtures.hpp"

TEST_F(epetraMultiVectorR9C4VecS9Fixture,
       MVecEpetraProductEigenVector){

  using namespace rompp;

  using mvec_t = core::MultiVector<Epetra_MultiVector>;
  STATIC_ASSERT_IS_CORE_MULTI_VECTOR_WRAPPER(mvec_t);
  mvec_t MV(*mv_);

  EXPECT_EQ( MV.globalNumVectors(), 4 );
  EXPECT_EQ( MV.localNumVectors(), 4 );
  EXPECT_EQ( MV.globalLength(), 9 );
  EXPECT_EQ( MV.localLength(), 3);
  for (int i=0; i<localSize_; i++)
    for (int j=0; j<MV.globalNumVectors(); j++)
      EXPECT_NEAR( 0.0, MV(i,j), 1e-12);

// MyPID    GID    
//     0     0      3.2     1.2     0      0
//     0     1      1.2      0      0      0
//     0     2       0       4      0      0
//     1     3       0       0      0      0
//     1     4       0       0      0      0
//     1     5       0       0      3      0
//     2     6       0       0      0      0
//     2     7       0       0      0      0
//     2     8       0       0      0      0
  
  if(rank_==0){
    MV(0,0) = 3.2;
    MV(1,0) = 1.2;
    MV(2,1) = 4;
    MV(0,1) = 1.2;}
  if(rank_==1){
    MV(2,2) = 3;}

  MV.data()->Print(std::cout);
  //----------

  using eigv_t = Eigen::Matrix<double,4,1>;
  eigv_t bn;
  using vec_t = core::Vector<eigv_t>;
  STATIC_ASSERT_IS_CORE_VECTOR_WRAPPER(vec_t);
  vec_t b(bn);
  b[0] = 1.;
  b[1] = 2.;
  b[2] = 3.;
  b[3] = 4.;

  auto res = core::ops::product(MV, b);
  res.data()->Print(std::cout);

  if (rank_==0){
    EXPECT_DOUBLE_EQ( res[0], 5.6);
    EXPECT_DOUBLE_EQ( res[1], 1.2);
    EXPECT_DOUBLE_EQ( res[2], 8.0);
  }

  if (rank_==1){
    EXPECT_DOUBLE_EQ( res[0], 0.0);
    EXPECT_DOUBLE_EQ( res[1], 0.0);
    EXPECT_DOUBLE_EQ( res[2], 9.0);
  }

  if (rank_==2){
    EXPECT_DOUBLE_EQ( res[0], 0.);
    EXPECT_DOUBLE_EQ( res[1], 0.);
    EXPECT_DOUBLE_EQ( res[2], 0.);
  }

}
