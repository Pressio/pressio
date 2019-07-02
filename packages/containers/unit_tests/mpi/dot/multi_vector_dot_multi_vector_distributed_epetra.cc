
#include "epetra_only_fixtures.hpp"

TEST(epetraMultiVectorR9C4VecS9Fixture,
       MVVecDotProduct){

  /*
    MultiVector A
    0 1 2 3
    1 0 1 2
    -------
    0 2 0 1
    0 2 2 1
    -------
    0 3 0 0
    4 0 1 0

    MultiVector B
    1 2
    0 0
    ----
    1 0
    1 0
    ----
    0 2
    0 1


    A^T B =  0 4
	     5 8
	     4 5
	     5 6
   */

  using namespace pressio;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  auto rank = comm.MyPID();
  EXPECT_EQ(comm.NumProc(),3);
  Epetra_Map rowMap(6, 0, comm);
  auto localSize = 2;

  using mvec_t = containers::MultiVector<Epetra_MultiVector>;
  mvec_t A(rowMap, 4);
  mvec_t B(rowMap, 2);
  static_assert(containers::meta::is_multi_vector_wrapper_epetra<mvec_t>::value,
		"");

  EXPECT_EQ( A.globalNumVectors(), 4 );
  EXPECT_EQ( B.localNumVectors(), 2 );
  EXPECT_EQ( A.globalLength(), 6 );
  EXPECT_EQ( B.globalLength(), 6 );

  for (int i=0; i<localSize; i++){
    for (int j=0; j<A.globalNumVectors(); j++)
      EXPECT_NEAR( 0.0, A(i,j), 1e-12);
    for (int j=0; j<B.globalNumVectors(); j++)
      EXPECT_NEAR( 0.0, B(i,j), 1e-12);
  }

  if(rank==0){
    A(0,0)=0.; A(0,1)=1.; A(0,2)=2.; A(0,3)=3.;
    A(1,0)=1.; A(1,1)=0.; A(1,2)=1.; A(1,3)=2.;
    B(0,0)=1.; B(0,1)=2.;
    B(1,0)=0.; B(1,1)=0.;
  }
  if(rank==1){
    A(0,0)=0.; A(0,1)=2.; A(0,2)=0.; A(0,3)=1.;
    A(1,0)=0.; A(1,1)=2.; A(1,2)=2.; A(1,3)=1.;
    B(0,0)=1.; B(0,1)=0.;
    B(1,0)=1.; B(1,1)=0.;
  }
  if(rank==2){
    A(0,0)=0.; A(0,1)=3.; A(0,2)=0.; A(0,3)=0.;
    A(1,0)=4.; A(1,1)=0.; A(1,2)=1.; A(1,3)=0.;
    B(0,0)=0.; B(0,1)=2.;
    B(1,0)=0.; B(1,1)=1.;
  }

  using eig_mat = Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic>;
  eig_mat TT(4,2);
  TT(0,0) = 0.0; TT(0,1) = 4.0;
  TT(1,0) = 5.0; TT(1,1) = 8.0;
  TT(2,0) = 4.0; TT(2,1) = 5.0;
  TT(3,0) = 5.0; TT(3,1) = 6.0;

  auto C = containers::ops::dot(A,B);
  for (auto i=0; i<4; i++){
    for (auto j=0; j<2; j++){
      EXPECT_NEAR( TT(i,j), C(i,j), 1e-12);
    }
  }
}
