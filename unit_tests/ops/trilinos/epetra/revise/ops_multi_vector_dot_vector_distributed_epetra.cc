
#include "epetra_only_fixtures.hpp"
#include "pressio_ops.hpp"

TEST_F(epetraMultiVectorR9C4VecS9Fixture,
       MVDotVec){

  using namespace pressio;

  assert(numProc_ == 3);
  using mvec_t = containers::MultiVector<Epetra_MultiVector>;
  using vec_t = containers::Vector<Epetra_Vector>;
  mvec_t MV(*mv_);
  vec_t b(*x_);

  EXPECT_EQ( MV.numVectors(), 4 );
  EXPECT_EQ( MV.numVectorsLocal(), 4 );
  EXPECT_EQ( MV.extent(0), 9 );
  EXPECT_EQ( MV.extentLocal(0), 3);

  for (int i=0; i<localSize_; i++)
    for (int j=0; j<MV.numVectors(); j++)
      EXPECT_NEAR( 0.0, MV(i,j), 1e-12);

  if(rank_==0){
    MV(0,0) = 3.2;
    MV(1,0) = 1.2;
    MV(2,1) = 4;
    MV(0,1) = 1.2;
  }

  if(rank_==1){
    MV(2,2) = 3;
  }

  b[0] = 1.0;
  b[1] = 1.0;
  b[2] = 1.0;

  // store into eigen dynamic vector wrapper
  //MV dot b = c2
  using natvec_t2 = Eigen::VectorXd;
  containers::Vector<natvec_t2> c2;
  c2.data()->resize(4);

  {  
  // constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
  // constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
  ::pressio::ops::product(::pressio::transpose(), 1., MV, b, 0., c2);
  EXPECT_EQ( c2.extent(0), 4);
  EXPECT_NEAR(c2[0], 4.4, 1e-12);
  EXPECT_NEAR(c2[1], 5.2, 1e-12);
  EXPECT_NEAR(c2[2], 3., 1e-12);
  EXPECT_NEAR(c2[3], 0., 1e-12);
  }

  // // store into teuchos serial dense vector wrapper
  // //MV dot b = c3
  // using natvec_t3 = Teuchos::SerialDenseVector<int, double>;
  // containers::Vector<natvec_t3> c3(4);
  // ops::dot(MV, b, c3);
  // EXPECT_EQ( c3.extent(0), 4);
  // EXPECT_NEAR(c3[0], 4.4, 1e-12);
  // EXPECT_NEAR(c3[1], 5.2, 1e-12);
  // EXPECT_NEAR(c3[2], 3., 1e-12);
  // EXPECT_NEAR(c3[3], 0., 1e-12);

}
