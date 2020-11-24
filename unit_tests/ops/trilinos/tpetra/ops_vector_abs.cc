
#include "tpetra_only_fixtures.hpp"
#include "pressio_ops.hpp"

TEST(TpetraOps, vectorAbs)
{
  {
    using tcomm = Teuchos::Comm<int>;
    using map_t = Tpetra::Map<>;
    using vec_t = Tpetra::Vector<>;

    Teuchos::RCP<const tcomm> comm =
      Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    auto numProc = comm->getSize();
    EXPECT_EQ(numProc,3);

    int numLocalEntries = 2;
    auto numGlobalEntries = numProc * numLocalEntries;
    Teuchos::RCP<const map_t> map =
      Teuchos::rcp(new map_t(numGlobalEntries, 0, comm));

    pressio::containers::Vector<vec_t> x(map);
    pressio::ops::fill(x, -3.1);
    pressio::containers::Vector<vec_t> y(map);
    pressio::ops::abs(y,x);

    auto yh = y.data()->getLocalViewHost();
    ASSERT_EQ( y.extent(0), 6 );
    for (int i=0; i<2; ++i) EXPECT_DOUBLE_EQ(yh(i,0), 3.1);
  }
}
