
#include "block_tpetra_only_fixtures.hpp"
#include "pressio_ops.hpp"

TEST(TpetraBlockOps, vectorAbs)
{
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using vec_t = Tpetra::BlockVector<>;

  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  auto numProc = comm->getSize();
  EXPECT_EQ(numProc,3);

  int numLocalEntries = 2;
  auto numGlobalEntries = numProc * numLocalEntries;
  Teuchos::RCP<const map_t> map =
    Teuchos::rcp(new map_t(numGlobalEntries, 0, comm));

  // 1. create vector x and fill with data
  pressio::containers::Vector<vec_t> x(*map, 2);
  pressio::ops::fill(x, -3.1);
  pressio::containers::Vector<vec_t> y(*map, 2);
  auto y_tpb = *y.data();
  pressio::ops::abs(y,x);

  ASSERT_EQ( y.extent(0), 6 );
  for (int i=0; i<numLocalEntries; ++i){
     auto lb = y_tpb.getLocalBlock(i);
     //std::cout << "*" << rank << " " << lb(0) << " " << lb(1) << std::endl;
     EXPECT_DOUBLE_EQ(lb[0], 3.1);
     EXPECT_DOUBLE_EQ(lb[1], 3.1);
  }
}
