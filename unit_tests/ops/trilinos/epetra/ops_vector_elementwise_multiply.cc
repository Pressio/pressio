
#include "epetra_only_fixtures.hpp"
#include "pressio_ops.hpp"

TEST(EpetraOps, elementwiseMultiplyVector)
{
  //using tcomm = Teuchos::Comm<int>;
  using map_t = Epetra_Map;
  using vec_t = Epetra_Vector;
  using ST = double;

  std::shared_ptr<Epetra_MpiComm> comm = std::make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
  auto rank    = comm->MyPID();
  auto numProc = comm->NumProc();
  EXPECT_EQ(numProc,3);

  int numLocalEntries = 2;
  auto numGlobalEntries = numProc * numLocalEntries;
  map_t map(numGlobalEntries, 0, *comm);

  auto shift = rank*numLocalEntries;

  // c = beta*c + alpha*x*b;

  // 1. create vector x and fill with data
  pressio::containers::Vector<vec_t> x(map);
  std::vector<ST> x0(numGlobalEntries);
  for (int i=0; i<6; ++i) x0[i]= (double) i;
  for (int i=0; i<numLocalEntries; ++i) x(i) = x0[shift+i];

  // 2. create vector b and fill with data
  pressio::containers::Vector<vec_t> b(map);
  std::vector<ST> b0(numGlobalEntries);
  for (int i=0; i<6; ++i) b0[i]= (double) i + 2.;
  for (int i=0; i<numLocalEntries; ++i) b(i) = b0[shift+i];

  // 3. compute ewise multiply
  pressio::containers::Vector<vec_t> c(map);
  pressio::ops::fill(c, 0.);
  pressio::ops::elementwise_multiply(1., x, b, 0., c);

  // 5. check correctness
  Eigen::VectorXd g(numGlobalEntries);
  g(0) = 0.;
  g(1) = 3.;
  g(2) = 8.;
  g(3) = 15.;
  g(4) = 24.;
  g(5) = 35.;

  ASSERT_EQ( c.extent(0), 6 );
  for (int i=0; i<numLocalEntries; ++i)
    EXPECT_DOUBLE_EQ(c(i), g(shift+i));
}
