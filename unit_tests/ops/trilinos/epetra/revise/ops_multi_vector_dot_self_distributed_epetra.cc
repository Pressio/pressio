
#include "epetra_only_fixtures.hpp"
#include "pressio_ops.hpp"

TEST(epetraMultiVector,
     MVVecDotSelf){

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

    result is symmetric
    A^T A =  17  0    5   2
	      0  18   6   7
	      5  6   10  10
	      2	 7   10  15
   */

  using namespace pressio;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  auto rank = comm.MyPID();
  assert(comm.NumProc() == 3);
  Epetra_Map rowMap(6, 0, comm);
  auto localSize = 2;

  using mvec_t = containers::MultiVector<Epetra_MultiVector>;
  mvec_t A(rowMap, 4);

  EXPECT_EQ( A.numVectors(), 4 );
  EXPECT_EQ( A.extent(0), 6 );

  for (int i=0; i<localSize; i++){
    for (int j=0; j<A.numVectors(); j++)
      EXPECT_NEAR( 0.0, A(i,j), 1e-12);
  }

  if(rank==0){
    A(0,0)=0.; A(0,1)=1.; A(0,2)=2.; A(0,3)=3.;
    A(1,0)=1.; A(1,1)=0.; A(1,2)=1.; A(1,3)=2.;
  }
  if(rank==1){
    A(0,0)=0.; A(0,1)=2.; A(0,2)=0.; A(0,3)=1.;
    A(1,0)=0.; A(1,1)=2.; A(1,2)=2.; A(1,3)=1.;
  }
  if(rank==2){
    A(0,0)=0.; A(0,1)=3.; A(0,2)=0.; A(0,3)=0.;
    A(1,0)=4.; A(1,1)=0.; A(1,2)=1.; A(1,3)=0.;
  }

  using eig_mat = Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic>;
  using eig_mat_w = containers::Matrix<eig_mat>;

  eig_mat TT(4,4);
  TT(0,0)=17.; TT(0,1)=0.; TT(0,2)=5.; TT(0,3)=2.;
  TT(1,0)=0.; TT(1,1)=18.; TT(1,2)=6.; TT(1,3)=7.;
  TT(2,0)=5.; TT(2,1)=6.; TT(2,2)=10.; TT(2,3)=10.;
  TT(3,0)=2.; TT(3,1)=7.; TT(3,2)=10.; TT(3,3)=15.;

  auto C = ops::product<eig_mat_w>(::pressio::transpose(), ::pressio::nontranspose(), 1., A);
  for (auto i=0; i<4; i++){
    for (auto j=0; j<4; j++){
      EXPECT_NEAR( TT(i,j), C(i,j), 1e-12);
    }
  }

  eig_mat_w C2(4,4);
  ops::product(::pressio::transpose(), ::pressio::nontranspose(), 1., A, 1., C2);
  for (auto i=0; i<4; i++){
    for (auto j=0; j<4; j++){
      EXPECT_NEAR( TT(i,j), C2(i,j), 1e-12);
    }
  }
}
