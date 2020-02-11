
#include "epetra_only_fixtures.hpp"

TEST_F(epetraMultiVectorR9C4VecS9Fixture,
       MVecEpetraProductEigenViewColVectorFromMV){

  using namespace pressio;

  using mvec_t = containers::MultiVector<Epetra_MultiVector>;
  STATIC_ASSERT_IS_CONTAINERS_MULTI_VECTOR_WRAPPER(mvec_t);
  mvec_t MV(*mv_);

  EXPECT_EQ( MV.numVectors(), 4 );
  EXPECT_EQ( MV.numVectorsLocal(), 4 );
  EXPECT_EQ( MV.extent(0), 9 );
  EXPECT_EQ( MV.extentLocal(0), 3);
  for (int i=0; i<localSize_; i++)
    for (int j=0; j<MV.numVectors(); j++)
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

  using eig_mat_t = Eigen::MatrixXd;
  eig_mat_t Bn(4,2);
  using mv_t = containers::MultiVector<eig_mat_t>;
  mv_t B(Bn);
  B(0, 0) = 1.; B(0, 1) = 2.;
  B(1, 0) = 2.; B(1, 1) = 4.;
  B(2, 0) = 3.; B(2, 1) = 6.;
  B(3, 0) = 4.; B(3, 1) = 8.;

  {
    // view the vector at j=0
    const auto colVec = pressio::containers::viewColumnVector(B, 0);
    EXPECT_DOUBLE_EQ( colVec[0], B(0,0) );
    EXPECT_DOUBLE_EQ( colVec[3], B(3,0) );
  }

  {
    // view the vector at j=1
    const auto colVec = pressio::containers::viewColumnVector(B, 1);
    EXPECT_DOUBLE_EQ( colVec[1], B(1,1) );
    EXPECT_DOUBLE_EQ( colVec[2], B(2,1) );
  }

  // do product
  {
    const auto colVec = pressio::containers::viewColumnVector(B, 0);
    auto res = containers::ops::product(MV, colVec);
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

    res.data()->PutScalar(0.);
    containers::ops::product(MV, colVec, res);
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


  // do product with second col
  {
    const auto colVec = pressio::containers::viewColumnVector(B, 1);
    auto res = containers::ops::product(MV, colVec);
    res.data()->Print(std::cout);
    if (rank_==0){
      EXPECT_DOUBLE_EQ( res[0], 11.2);
      EXPECT_DOUBLE_EQ( res[1], 2.4);
      EXPECT_DOUBLE_EQ( res[2], 16.0);
    }
    if (rank_==1){
      EXPECT_DOUBLE_EQ( res[0], 0.0);
      EXPECT_DOUBLE_EQ( res[1], 0.0);
      EXPECT_DOUBLE_EQ( res[2], 18.0);
    }
    if (rank_==2){
      EXPECT_DOUBLE_EQ( res[0], 0.);
      EXPECT_DOUBLE_EQ( res[1], 0.);
      EXPECT_DOUBLE_EQ( res[2], 0.);
    }

    res.data()->PutScalar(0.);
    containers::ops::product(MV, colVec, res);
    if (rank_==0){
      EXPECT_DOUBLE_EQ( res[0], 11.2);
      EXPECT_DOUBLE_EQ( res[1], 2.4);
      EXPECT_DOUBLE_EQ( res[2], 16.0);
    }
    if (rank_==1){
      EXPECT_DOUBLE_EQ( res[0], 0.0);
      EXPECT_DOUBLE_EQ( res[1], 0.0);
      EXPECT_DOUBLE_EQ( res[2], 18.0);
    }
    if (rank_==2){
      EXPECT_DOUBLE_EQ( res[0], 0.);
      EXPECT_DOUBLE_EQ( res[1], 0.);
      EXPECT_DOUBLE_EQ( res[2], 0.);
    }
  }

}
