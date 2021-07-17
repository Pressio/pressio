
#include "block_tpetra_only_fixtures.hpp"
#include "pressio_ops.hpp"

TEST(TpetraBlockOps, elementwiseMultiplyVector)
{
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using vec_t = Tpetra::BlockVector<>;
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
  pressio::containers::Vector<vec_t> x(*map, 2);
  auto & x_tpb = *x.data();
  std::vector<ST> x0(numGlobalEntries);
  for (int i=0; i<6; ++i) x0[i]= (double) i;
  for (int i=0; i<numLocalEntries; ++i){
    ST vals[3];
    vals[0] = x0[shift+i];
    vals[1] = x0[shift+i];
    x_tpb.replaceLocalValues(i, vals);
    //auto lb0 = x_tpb.getLocalBlock(i);
    //std::cout << rank << " " << i << " " << lb0(0) << " " << lb0(1) << std::endl;
  }

  // 2. create vector b and fill with data
  pressio::containers::Vector<vec_t> b(*map, 2);
  auto & b_tpb = *b.data();
  std::vector<ST> b0(numGlobalEntries);
  for (int i=0; i<6; ++i) b0[i]= (double) i + 2.;
  for (int i=0; i<numLocalEntries; ++i){
    ST vals[3];
    vals[0] = b0[shift+i];
    vals[1] = b0[shift+i];
    b_tpb.replaceLocalValues(i, vals);
    //auto lb0 = b_tpb.getLocalBlock(i);
  }

  // 3. compute ewise multiply
  pressio::containers::Vector<vec_t> c(*map, 2);
  pressio::ops::fill(c, 0.);
  pressio::ops::elementwise_multiply(1., x, b, 0., c);
  // auto ch = c.data()->getLocalViewHost();

  // 5. check correctness
  Eigen::VectorXd g(numGlobalEntries);
  g(0) = 0.;
  g(1) = 3.;
  g(2) = 8.;
  g(3) = 15.;
  g(4) = 24.;
  g(5) = 35.;
  ASSERT_EQ( c.extent(0), 6 );
  for (int i=0; i<numLocalEntries; ++i){
    auto lb = c.data()->getLocalBlock(i);
     EXPECT_DOUBLE_EQ(lb[0], g(shift+i));
     EXPECT_DOUBLE_EQ(lb[1], g(shift+i));
  }
}
