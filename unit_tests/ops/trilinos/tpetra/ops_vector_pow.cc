
#include "tpetra_only_fixtures.hpp"
#include "pressio_ops.hpp"

TEST(TpetraOps, vectorPow)
{
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using vec_t = Tpetra::Vector<>;
  using ST = typename vec_t::scalar_type;

  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  auto rank    = comm->getRank();
  auto numProc = comm->getSize();
  EXPECT_EQ(numProc,3);

  int numLocalEntries = 2;
  auto numGlobalEntries = numProc * numLocalEntries;
  Teuchos::RCP<const map_t> map =
    Teuchos::rcp(new map_t(numGlobalEntries, 0, comm));

  auto shift = rank*numLocalEntries;

  // 1. create vector x and fill with data
  pressio::containers::Vector<vec_t> x(map);
  auto xh = x.data()->getLocalViewHost();
  std::vector<ST> x0(numGlobalEntries);
  for (int i=0; i<6; ++i) x0[i]= (double) i;
  for (int i=0; i<numLocalEntries; ++i) xh(i,0) = x0[shift+i];

  // 2. compute pow
  ::pressio::ops::pow(x, 2.);

  // 5. check correctness
  Eigen::VectorXd g(numGlobalEntries);
  g(0) = 0.;
  g(1) = 1.;
  g(2) = 4.;
  g(3) = 9.;
  g(4) = 16.;
  g(5) = 25.;
  ASSERT_EQ( x.extent(0), 6 );
  for (int i=0; i<numLocalEntries; ++i)
    EXPECT_DOUBLE_EQ(xh(i,0), g(shift+i));
}
