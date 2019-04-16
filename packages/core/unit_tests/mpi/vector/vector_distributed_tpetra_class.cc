
#include "tpetra_only_fixtures.hpp"

using native_t = typename tpetraVectorGlobSize15Fixture::vec_t;


TEST_F(tpetraVectorGlobSize15Fixture, Constructor){
  using namespace rompp;

  std::cout << rank_ << std::endl;
  using myvec_t = core::Vector<native_t>;
  myvec_t a( *x_ );
  myvec_t a3( a );
}

TEST_F(tpetraVectorGlobSize15Fixture,
       localSize){
  using namespace rompp;
  using myvec_t = core::Vector<native_t>;
  myvec_t v1( *x_ );
  EXPECT_EQ( v1.localSize(), 5);
}

TEST_F(tpetraVectorGlobSize15Fixture,
       globalSize){
  using namespace rompp;
  using myvec_t = core::Vector<native_t>;
  myvec_t v1( *x_ );
  EXPECT_EQ( v1.globalSize(), 15);
}

TEST_F(tpetraVectorGlobSize15Fixture,
       isGloballyDist){
  using namespace rompp;
  using myvec_t = core::Vector<native_t>;
  myvec_t v1( *x_ );
  EXPECT_TRUE( v1.isDistributedGlobally() );
}

TEST_F(tpetraVectorGlobSize15Fixture,
       SetScalar){
  using namespace rompp;
  using sc_t = typename tpetraVectorGlobSize15Fixture::ST;

  using myvec_t = core::Vector<native_t>;
  myvec_t v1( *x_ );
  v1.putScalar(43.3);

  Teuchos::ArrayRCP<const sc_t> dd = v1.data()->getData();
  for (int i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( dd[i], 43.3 );
  }
}

TEST_F(tpetraVectorGlobSize15Fixture,
       CompoundAssignAdd_deep_copy_cstr){
  using namespace rompp;
  using sc_t = typename tpetraVectorGlobSize15Fixture::ST;

  using myvec_t = core::Vector<native_t>;
  myvec_t v1( *x_ );
  v1.putScalar(1.2);
  myvec_t v2( *x_ );
  v2.putScalar(3.3);
  v1 += v2;
  Teuchos::ArrayRCP<const sc_t> dd = v1.data()->getData();
  for (int i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( dd[i], 4.5 );
  }
}

TEST_F(tpetraVectorGlobSize15Fixture,
       CompoundAssignAdd_mapConstr){
  using namespace rompp;
  using sc_t = typename tpetraVectorGlobSize15Fixture::ST;

  using myvec_t = core::Vector<native_t>;
  myvec_t v1( contigMap_ );
  v1.putScalar(1.2);
  myvec_t v2( contigMap_ );
  v2.putScalar(3.3);
  v1 += v2;
  Teuchos::ArrayRCP<const sc_t> dd = v1.data()->getData();
  for (int i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( dd[i], 4.5 );
  }
}

TEST_F(tpetraVectorGlobSize15Fixture,
       CompoundAssignSubtract_deep_copy_cstr){
  using namespace rompp;
  using sc_t = typename tpetraVectorGlobSize15Fixture::ST;

  using myvec_t = core::Vector<native_t>;
  myvec_t v1( *x_ );
  v1.putScalar(1.2);
  myvec_t v2( *x_ );
  v2.putScalar(3.3);
  v1 -= v2;
  Teuchos::ArrayRCP<const sc_t> dd = v1.data()->getData();
  for (int i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( dd[i], -2.1 );
  }
}

TEST_F(tpetraVectorGlobSize15Fixture,
       CompoundAssignSubtract_mapConstr){
  using namespace rompp;
  using sc_t = typename tpetraVectorGlobSize15Fixture::ST;

  using myvec_t = core::Vector<native_t>;
  myvec_t v1( contigMap_ );
  v1.putScalar(1.2);
  myvec_t v2( contigMap_ );
  v2.putScalar(3.3);
  v1 -= v2;
  Teuchos::ArrayRCP<const sc_t> dd = v1.data()->getData();
  for (int i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( dd[i], -2.1 );
  }
}

TEST_F(tpetraVectorGlobSize15Fixture,
       SetZero){
  using namespace rompp;
  using sc_t = typename tpetraVectorGlobSize15Fixture::ST;

  using myvec_t = core::Vector<native_t>;
  myvec_t v1( *x_ );
  v1.setZero();

  Teuchos::ArrayRCP<const sc_t> dd = v1.data()->getData();
  for (int i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( dd[i], 0.0 );
  }
}

TEST_F(tpetraVectorGlobSize15Fixture,
       QueryWrappedData){
  using namespace rompp;

  using myvec_t = core::Vector<native_t>;
  myvec_t v1( *x_ );
  ::testing::StaticAssertTypeEq<decltype(v1.data()),
  				native_t * >();
  const myvec_t v2( *x_ );
  ::testing::StaticAssertTypeEq< decltype(v2.data()),
  				 const native_t * >();
}

TEST_F(tpetraVectorGlobSize15Fixture,
       empty){
  using namespace rompp;
  using myvec_t = core::Vector<native_t>;
  myvec_t v1( *x_ );
  EXPECT_FALSE(v1.empty());
}

TEST_F(tpetraVectorGlobSize15Fixture,
       getMap){
  using namespace rompp;
  using myvec_t = core::Vector<native_t>;
  myvec_t v1( *x_ );
  auto const & mapO = v1.getDataMap();

  ::testing::StaticAssertTypeEq<decltype(mapO),
  				const typename tpetraVectorGlobSize15Fixture::map_t & >();
  EXPECT_TRUE(mapO.isContiguous());
}
