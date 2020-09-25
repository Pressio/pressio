
#include "epetra_only_fixtures.hpp"

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
  for (int i=0; i<v1.extentLocal(0); i++){
    v1[i] = 11.2;
  }
  for (int i=0; i<v1.extentLocal(0); i++){
    EXPECT_DOUBLE_EQ( v1[i], (*x_)[i] );
  }
  v1[3] = 56.;
  EXPECT_DOUBLE_EQ( v1[3], 56.0);
}

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

  // changing b does not change a
  b.data()->PutScalar(2.);
  b.data()->Norm1(&res); ASSERT_EQ( res, 30.);
  a.Norm1(&res); ASSERT_EQ( res, 15.);

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
  // check that b has correct values
  double res; b.data()->Norm1(&res);
  ASSERT_EQ( res, 15. );

  // changing b does not change a
  b.data()->PutScalar(2.);
  b.data()->Norm1(&res); ASSERT_EQ( res, 30.);
  a.data()->Norm1(&res); ASSERT_EQ( res, 15.);

  // check underlying ptr is different than a's ptr and both are valid
  ASSERT_TRUE( a.data()->Values() != b.data()->Values() );
  ASSERT_TRUE( a.data()->Values() != nullptr );
  ASSERT_TRUE( b.data()->Values() != nullptr );
}

TEST(containers_vector_epetra, MoveConstr)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  Epetra_Map map(15, 0, comm);

  w_t a(map);
  a.data()->PutScalar(1.);

  w_t b(std::move(a));
  // check that b has correct values
  double res; b.data()->Norm1(&res);
  ASSERT_EQ( res, 15. );

  // epetra does not have move semnatics so should behave like copy
  ASSERT_TRUE( a.data()->Values() != b.data()->Values() );
  ASSERT_TRUE( a.data()->Values() != nullptr );
  ASSERT_TRUE( b.data()->Values() != nullptr );
}

TEST(containers_vector_epetra, MoveAssign)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  Epetra_Map map(15, 0, comm);

  w_t b(map);
  {
    w_t a(map);
    a.data()->PutScalar(1.);
    b = std::move(a);
    // epetra does not have move semnatics so should behave like copy
    ASSERT_TRUE( a.data()->Values() != b.data()->Values() );
    ASSERT_TRUE( a.data()->Values() != nullptr );
    ASSERT_TRUE( b.data()->Values() != nullptr );
  }

  // check that b has correct values
  double res; b.data()->Norm1(&res);
  ASSERT_EQ( res, 15. );
}
