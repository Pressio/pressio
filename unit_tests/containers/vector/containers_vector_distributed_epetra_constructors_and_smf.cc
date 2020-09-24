
#include <gtest/gtest.h>
#include "Epetra_MpiComm.h"
#include "pressio_containers.hpp"

using vec_t = Epetra_Vector;
using w_t = pressio::containers::Vector<vec_t>;

TEST(containers_vector_epetra, Constructor1)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  Epetra_Map map(15, 0, comm);
  w_t b(map);
  ASSERT_TRUE( b.data()->Values() != nullptr );
  ASSERT_EQ( b.extent(0),     15);
  ASSERT_EQ( b.extentLocal(0), 5);
}

TEST(containers_vector_epetra, Constructor2)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  Epetra_Map map(15, 0, comm);
  Epetra_Vector a(map);
  a.PutScalar(1.);

  w_t b(a);
  double res; b.data()->Norm1(&res);
  ASSERT_EQ( res, 15. );
  ASSERT_TRUE( b.data()->Values() != nullptr );
  ASSERT_TRUE( b.data()->Values() != a.Values() );
  ASSERT_TRUE( b.data()->Values() != nullptr );
}

TEST(containers_vector_epetra, CopyConstructor)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  Epetra_Map map(15, 0, comm);

  w_t a(map);
  a.data()->PutScalar(1.);

  w_t b(a);
  double res; b.data()->Norm1(&res);
  ASSERT_EQ( res, 15. );
  ASSERT_TRUE( a.data()->Values() != b.data()->Values() );
  ASSERT_TRUE( a.data()->Values() != nullptr );
  ASSERT_TRUE( b.data()->Values() != nullptr );
}

// TEST(containers_vector_epetra, CopyAssign)
// {
//   Epetra_MpiComm comm(MPI_COMM_WORLD);
//   Epetra_Map map(15, 0, comm);

//   w_t a(map);
//   a.data()->PutScalar(1.);

//   w_t b(map);
//   b = a;
//   double res; b.data()->Norm1(&res);
//   ASSERT_EQ( res, 15. );
//   ASSERT_TRUE( a.data()->Values() != b.data()->Values() );
//   ASSERT_TRUE( a.data()->Values() != nullptr );
//   ASSERT_TRUE( b.data()->Values() != nullptr );
// }

TEST(containers_vector_epetra, MoveConstr)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  Epetra_Map map(15, 0, comm);

  w_t a(map);
  a.data()->PutScalar(1.);

  w_t b(std::move(a));
  double res; b.data()->Norm1(&res);
  ASSERT_EQ( res, 15. );
  // epetra does not have move semantics
  ASSERT_TRUE( a.data()->Values() != b.data()->Values() );
  ASSERT_TRUE( a.data()->Values() != nullptr );
  ASSERT_TRUE( b.data()->Values() != nullptr );
}

TEST(containers_vector_epetra, MoveAssign)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  Epetra_Map map(15, 0, comm);

  w_t a(map);
  a.data()->PutScalar(1.);

  w_t b(std::move(a));
  double res; b.data()->Norm1(&res);
  ASSERT_EQ( res, 15. );
  // epetra does not have move semantics
  ASSERT_TRUE( a.data()->Values() != b.data()->Values() );
  ASSERT_TRUE( a.data()->Values() != nullptr );
  ASSERT_TRUE( b.data()->Values() != nullptr );
}
