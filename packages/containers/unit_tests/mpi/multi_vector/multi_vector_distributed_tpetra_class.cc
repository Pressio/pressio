
#include "tpetra_only_fixtures.hpp"

TEST_F(tpetraMultiVectorGlobSize15Fixture, Constructor){
  using namespace rompp;

  std::cout << rank_ << std::endl;
  using mvec_t = containers::MultiVector<typename tpetraMultiVectorGlobSize15Fixture::mvec_t>;
  mvec_t a( *x_ );
  mvec_t a3( a );
}


TEST_F(tpetraMultiVectorGlobSize15Fixture,
       localSize){
  using namespace rompp;
  using mymvec_t = containers::MultiVector<typename tpetraMultiVectorGlobSize15Fixture::mvec_t>;
  mymvec_t v1( *x_ );
  EXPECT_EQ( v1.localLength(), 5);
}

// TEST_F(tpetraMultiVectorGlobSize15Fixture,
//        print){
//   using namespace rompp;
//   using mymvec_t = containers::MultiVector<typename tpetraMultiVectorGlobSize15Fixture::mvec_t>;
//   mymvec_t v1( *x_ );
//   v1.print();
// }


TEST_F(tpetraMultiVectorGlobSize15Fixture,
       comm){
  using namespace rompp;
  using native_t = typename tpetraMultiVectorGlobSize15Fixture::mvec_t;
  using mymvec_t = containers::MultiVector<native_t>;
  mymvec_t v1( *x_ );

  auto commRCP = v1.comm();
  ::testing::StaticAssertTypeEq<decltype(commRCP),
      Teuchos::RCP< const Teuchos::Comm<int>> >();
}


TEST_F(tpetraMultiVectorGlobSize15Fixture,
       globalSize){
  using namespace rompp;
  using mymvec_t = containers::MultiVector<typename tpetraMultiVectorGlobSize15Fixture::mvec_t>;
  mymvec_t v1( *x_ );
  EXPECT_EQ( v1.globalLength(), 15);
}


TEST_F(tpetraMultiVectorGlobSize15Fixture,
       isGloballyDist){
  using namespace rompp;
  using mymvec_t = containers::MultiVector<typename tpetraMultiVectorGlobSize15Fixture::mvec_t>;
  mymvec_t v1( *x_ );
  EXPECT_TRUE( v1.isDistributedGlobally() );
}


TEST_F(tpetraMultiVectorGlobSize15Fixture,
       QueryWrappedData){
  using namespace rompp;
  using nvec_t = typename tpetraMultiVectorGlobSize15Fixture::mvec_t;

  using mymvec_t = containers::MultiVector<nvec_t>;
  mymvec_t v1( *x_ );
  ::testing::StaticAssertTypeEq<decltype(v1.data()),
  				nvec_t * >();
  const mymvec_t v2( *x_ );
  ::testing::StaticAssertTypeEq< decltype(v2.data()),
  				 const nvec_t * >();
}


TEST_F(tpetraMultiVectorGlobSize15Fixture,
       empty){
  using namespace rompp;
  using nvec_t = typename tpetraMultiVectorGlobSize15Fixture::mvec_t;
  using mymvec_t = containers::MultiVector<nvec_t>;
  mymvec_t v1( *x_ );
  EXPECT_FALSE(v1.empty());
}


// TEST_F(tpetraMultiVectorGlobSize15Fixture,
//        getMap){
//   using namespace rompp;
//   using nvec_t = typename tpetraMultiVectorGlobSize15Fixture::mvec_t;
//   using mymvec_t = containers::MultiVector<nvec_t>;
//   mymvec_t v1( *x_ );
//   auto const & mapO = v1.getDataMap();

//   ::testing::StaticAssertTypeEq<decltype(mapO),
//   				const typename tpetraMultiVectorGlobSize15Fixture::map_t & >();
//   EXPECT_TRUE(mapO.isContiguous());
// }


TEST_F(tpetraMultiVectorGlobSize15Fixture,
       localNumVecs){
  using namespace rompp;
  using nvec_t = typename tpetraMultiVectorGlobSize15Fixture::mvec_t;
  using mymvec_t = containers::MultiVector<nvec_t>;
  mymvec_t v1( *x_ );
  EXPECT_EQ(v1.localNumVectors(), 4);
}


TEST_F(tpetraMultiVectorGlobSize15Fixture,
       globalNumVecs){
  using namespace rompp;
  using nvec_t = typename tpetraMultiVectorGlobSize15Fixture::mvec_t;
  using mymvec_t = containers::MultiVector<nvec_t>;
  mymvec_t v1( *x_ );
  EXPECT_EQ(v1.globalNumVectors(), 4);
}


TEST_F(tpetraMultiVectorGlobSize15Fixture,
       SetZero){
  using namespace rompp;
  using sc_t = typename tpetraMultiVectorGlobSize15Fixture::ST;

  using mymvec_t = containers::MultiVector<typename tpetraMultiVectorGlobSize15Fixture::mvec_t>;
  mymvec_t v1( *x_ );
  v1.setZero();

  for (int k=0; k<4; k++){
    Teuchos::ArrayRCP<const sc_t> dd = v1.data()->getData(k);
    for (int i=0; i<v1.localLength(); i++){
      EXPECT_DOUBLE_EQ( dd[i], 0.0 );
    }
  }
}


TEST_F(tpetraMultiVectorGlobSize15Fixture,
       getMap){
  using namespace rompp;
  using nvec_t = typename tpetraMultiVectorGlobSize15Fixture::mvec_t;
  using myvec_t = containers::MultiVector<nvec_t>;
  myvec_t v1( *x_ );
  auto const & mapO = v1.getDataMap();
  ::testing::StaticAssertTypeEq<decltype(mapO),
  				const typename tpetraMultiVectorGlobSize15Fixture::map_t & >();
  EXPECT_TRUE(mapO.isContiguous());

  auto mapO1 = v1.getRCPDataMap();
  ::testing::StaticAssertTypeEq<decltype(mapO1),
  	Teuchos::RCP< const typename tpetraMultiVectorGlobSize15Fixture::map_t>>();
  EXPECT_TRUE(mapO1->isContiguous());
}


TEST_F(tpetraMultiVectorGlobSize15Fixture,
       constrcutTVectorFromMVMap){
  using namespace rompp;
  using nvec_t = typename tpetraMultiVectorGlobSize15Fixture::mvec_t;
  using myvec_t = containers::MultiVector<nvec_t>;
  myvec_t v1( *x_ );

  auto mapO1 = v1.getRCPDataMap();

  using sc_t = typename tpetraMultiVectorGlobSize15Fixture::ST;
  using LO_t = typename tpetraMultiVectorGlobSize15Fixture::LO;
  using GO_t = typename tpetraMultiVectorGlobSize15Fixture::GO;
  using NO_t = typename tpetraMultiVectorGlobSize15Fixture::NT;

  Tpetra::Vector<sc_t,LO_t,GO_t,NO_t> vec(mapO1);
  // ::testing::StaticAssertTypeEq<decltype(mapO1),
  // 	Teuchos::RCP< const typename tpetraMultiVectorGlobSize15Fixture::map_t>>();
  // EXPECT_TRUE(mapO1->isContiguous());
}
