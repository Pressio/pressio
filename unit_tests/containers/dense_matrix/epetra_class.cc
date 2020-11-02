
#include "epetra_only_fixtures.hpp"

using n_t = Epetra_MultiVector;
using w_t = pressio::containers::DenseMatrix<n_t>;

TEST(containers_mv_epetra, Constructor1)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  Epetra_Map map(15, 0, comm);
  w_t b(map, 3);
  ASSERT_TRUE( b.data()->Values() != nullptr );
  ASSERT_EQ( b.extent(0),     15);
  ASSERT_EQ( b.extent(1),     3);
  ASSERT_EQ( b.extentLocal(0), 5);
}

TEST(containers_mv_epetra, Constructor2)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  Epetra_Map map(15, 0, comm);
  n_t a(map, 3);
  a.PutScalar(1.);

  w_t b(a);
  double res[3]; b.data()->Norm1(res);
  ASSERT_EQ( res[0], 15. );

  // changing b does not change a
  b.data()->PutScalar(2.);
  b.data()->Norm1(res); 
  ASSERT_EQ( res[0], 30.);
  ASSERT_EQ( res[1], 30.);
  ASSERT_EQ( res[2], 30.);
  a.Norm1(res);         
  ASSERT_EQ( res[0], 15.);
  ASSERT_EQ( res[1], 15.);
  ASSERT_EQ( res[2], 15.);

  ASSERT_TRUE( b.data()->Values() != nullptr );
  ASSERT_TRUE( b.data()->Values() != a.Values() );
  ASSERT_TRUE( b.data()->Values() != nullptr );
}

TEST(containers_mv_epetra, CopyConstructor)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  Epetra_Map map(15, 0, comm);

  w_t a(map, 3);
  a.data()->PutScalar(1.);

  w_t b(a);
  // check that b has correct values
  double res[3]; b.data()->Norm1(res);
  ASSERT_EQ( res[0], 15.);
  ASSERT_EQ( res[1], 15.);
  ASSERT_EQ( res[2], 15.);

  // changing b does not change a
  b.data()->PutScalar(2.);
  b.data()->Norm1(res); 
  ASSERT_EQ( res[0], 30.);
  ASSERT_EQ( res[1], 30.);
  ASSERT_EQ( res[2], 30.);
  a.data()->Norm1(res); 
  ASSERT_EQ( res[0], 15.);
  ASSERT_EQ( res[1], 15.);
  ASSERT_EQ( res[2], 15.);

  // check underlying ptr is different than a's ptr and both are valid
  ASSERT_TRUE( a.data()->Values() != b.data()->Values() );
  ASSERT_TRUE( a.data()->Values() != nullptr );
  ASSERT_TRUE( b.data()->Values() != nullptr );
}

TEST(containers_mv_epetra, MoveConstr)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  Epetra_Map map(15, 0, comm);

  w_t a(map, 3);
  a.data()->PutScalar(1.);

  w_t b(std::move(a));
  // check that b has correct values
  double res[3]; 
  b.data()->Norm1(res);
  ASSERT_EQ( res[0], 15.);
  ASSERT_EQ( res[1], 15.);
  ASSERT_EQ( res[2], 15.);

  // epetra does not have move semnatics so should behave like copy
  ASSERT_TRUE( a.data()->Values() != b.data()->Values() );
  ASSERT_TRUE( a.data()->Values() != nullptr );
  ASSERT_TRUE( b.data()->Values() != nullptr );
}

TEST(containers_mv_epetra, MoveAssign)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  Epetra_Map map(15, 0, comm);

  w_t b(map, 3);
  {
    w_t a(map, 3);
    a.data()->PutScalar(1.);
    b = std::move(a);
    // epetra does not have move semnatics so should behave like copy
    ASSERT_TRUE( a.data()->Values() != b.data()->Values() );
    ASSERT_TRUE( a.data()->Values() != nullptr );
    ASSERT_TRUE( b.data()->Values() != nullptr );
  }

  // check that b has correct values
  double res[3]; b.data()->Norm1(res);
  ASSERT_EQ( res[0], 15.);
  ASSERT_EQ( res[1], 15.);
  ASSERT_EQ( res[2], 15.);  
}
