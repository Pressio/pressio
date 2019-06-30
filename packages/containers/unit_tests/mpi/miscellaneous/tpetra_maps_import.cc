
#include "tpetra_only_fixtures.hpp"

TEST_F(tpetraMultiVectorGlobSize9Fixture,
       nonContigMapsMimicCollocation){
  using namespace rompp;

  // aliases
  using nat_mv_t = typename tpetraMultiVectorGlobSize9Fixture::mvec_t;
  using mvec_t = containers::MultiVector<nat_mv_t>;
  using device_t = typename containers::details::traits<mvec_t>::device_t;
  using node_t = typename nat_mv_t::node_type;

  // construct multivector wrapper
  mvec_t MV( *x_ );
  MV.setZero();
  // get trilinos tpetra multivector object
  auto trilD = MV.data();
  trilD->sync<Kokkos::HostSpace>();

  auto v2d = trilD->getLocalView<Kokkos::HostSpace>();
  //we are going to change the host view
  trilD->modify<Kokkos::HostSpace>();

  // MyPID    GID
  //     0     0      3.2     1.2     0      0
  //     0     1      1.2      0      0      0
  //     0     2       0       4      0      0

  //     1     3       0       0      0      0
  //     1     4       0       5     -1      0
  //     1     5       0       0      3      0

  //     2     6       0       0      0      0
  //     2     7       0       0      0      44
  //     2     8       0       0      0      0

  if(rank_==0){
    v2d(0,0) = 3.2; v2d(0,1) = 1.2;
    v2d(1,0) = 1.2;
    v2d(2,1) = 4;
  }
  if(rank_==1){
    v2d(1,1) = 5; v2d(1,2) = -1;
    v2d(2,2) = 3;
  }
  if(rank_==2){
    v2d(1,3) = 44;
  }
  // sync from host to device
  trilD->sync<device_t> ();

  EXPECT_EQ( MV.globalNumVectors(), 4 );
  EXPECT_EQ( MV.localNumVectors(), 4 );
  EXPECT_EQ( MV.globalLength(), 9 );
  EXPECT_EQ( MV.localLength(), 3);

  //-------------------------------------
  //-------------------------------------

  // create non contig map
  std::vector<GO> myGID;
  if(rank_==0)
    myGID = {1};
  if(rank_==1)
    myGID = {4,5};
  if(rank_==2)
    myGID = {7};

  Teuchos::ArrayView< const GO > elL(myGID);
  //create a non contiguous map
  Teuchos::RCP<const Tpetra::Map<LO,GO>> nonCmap =
    Tpetra::createNonContigMap<LO,GO>(elL, comm_);
  // create the multivector for the non-contig map
  std::shared_ptr<nat_mv_t> A2 = std::make_shared<nat_mv_t>(nonCmap, numVecs_);

  using import_t = Tpetra::Import<LO, GO, node_t>;
  //                 (source,   target)
  import_t importer(contigMap_, nonCmap);
  A2->doImport(*MV.data(), importer, Tpetra::INSERT);
  // print
  Tpetra::MatrixMarket::Writer<nat_mv_t>::writeDense(std::cout, *A2, "void", "dfdfd");

  // check data
  A2->sync<Kokkos::HostSpace>();
  auto vv = A2->getLocalView<Kokkos::HostSpace>();

  if(rank_==0){
    EXPECT_DOUBLE_EQ(vv(0,0), 1.2);
  }
  if(rank_==1){
    EXPECT_DOUBLE_EQ(vv(0,1), 5);
    EXPECT_DOUBLE_EQ(vv(0,2), -1);
  }
  if(rank_==2){
    EXPECT_DOUBLE_EQ(vv(0,3), 44);
  }
}
