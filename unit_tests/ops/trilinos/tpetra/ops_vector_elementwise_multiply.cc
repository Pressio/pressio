
#include "tpetra_only_fixtures.hpp"
#include "pressio_ops.hpp"

TEST(TpetraOps, elementwiseMultiplyVector)
{
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using vec_t = Tpetra::Vector<>;
  using ST = typename vec_t::scalar_type;

  Teuchos::RCP<const tcomm> comm = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  auto rank    = comm->getRank();
  auto numProc = comm->getSize();
  EXPECT_EQ(numProc,3);

  int numLocalEntries = 2;
  auto numGlobalEntries = numProc * numLocalEntries;
  Teuchos::RCP<const map_t> map = Teuchos::rcp(new map_t(numGlobalEntries, 0, comm));

  auto shift = rank*numLocalEntries;

  // c = beta*c + alpha*x*b;

  // 1. create vector x and fill with data
  pressio::containers::Vector<vec_t> x(map);
  auto xh = x.data()->getLocalViewHost();
  std::vector<ST> x0(numGlobalEntries);
  for (int i=0; i<6; ++i) x0[i]= (double) i;
  for (int i=0; i<numLocalEntries; ++i) xh(i,0) = x0[shift+i];

  // 2. create vector b and fill with data
  pressio::containers::Vector<vec_t> b(map);
  auto bh = b.data()->getLocalViewHost();
  std::vector<ST> b0(numGlobalEntries);
  for (int i=0; i<6; ++i) b0[i]= (double) i + 2.;
  for (int i=0; i<numLocalEntries; ++i) bh(i,0) = b0[shift+i];

  // 3. compute ewise multiply
  pressio::containers::Vector<vec_t> c(map);
  pressio::ops::fill(c, 0.);
  pressio::ops::elementwise_multiply(1., x, b, 0., c);
  auto ch = c.data()->getLocalViewHost();

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
    EXPECT_DOUBLE_EQ(ch(i,0), g(shift+i));
}
