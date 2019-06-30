
#include "block_tpetra_only_fixtures.hpp"

using fix_t = tpetraBlockVectorGlobSize15BlockSize5Fixture;
using native_t = typename fix_t::vec_t;
using myvec_t = rompp::containers::Vector<native_t>;

TEST_F(tpetraBlockVectorGlobSize15BlockSize5Fixture,
       QueryWrappedData){
  using namespace rompp;

  using myvec_t = containers::Vector<native_t>;
  myvec_t v1( *x_ );
  ::testing::StaticAssertTypeEq<decltype(v1.data()),
  				native_t * >();
  const myvec_t v2( *x_ );
  ::testing::StaticAssertTypeEq< decltype(v2.data()),
  				 const native_t * >();
}

TEST_F(tpetraBlockVectorGlobSize15BlockSize5Fixture,
       getMap){
  using namespace rompp;
  using myvec_t = containers::Vector<native_t>;
  myvec_t v1( *x_ );
  auto const & mapO = v1.getDataMap();
  ::testing::StaticAssertTypeEq
      <decltype(mapO),
       const typename fix_t::map_t & >();
  EXPECT_TRUE(mapO.isContiguous());
}

TEST_F(tpetraBlockVectorGlobSize15BlockSize5Fixture,
       ConstructorFromNative){
  using namespace rompp;
  STATIC_ASSERT_IS_VECTOR_TPETRA_BLOCK(native_t);
  static_assert
    ( containers::meta::is_vector_wrapper_tpetra_block<myvec_t>::value, "");
  using sc_t = typename fix_t::ST;

  x_->putScalar(1.22);
  myvec_t a( *x_ );

  // a is wrapper of a block vector, so get a tpetra vector
  auto a_tp = a.data()->getVectorView();
  // now that we have a regular tpetra vector, we can check data
  Teuchos::ArrayRCP<const sc_t> dd = a_tp.getData();
  for (int i=0; i<a.localSize(); i++){
    EXPECT_DOUBLE_EQ( dd[i], 1.22 );
  }
}

TEST_F(tpetraBlockVectorGlobSize15BlockSize5Fixture,
       CopyConstructor){
  using namespace rompp;
  STATIC_ASSERT_IS_VECTOR_TPETRA_BLOCK(native_t);
  static_assert
    ( containers::meta::is_vector_wrapper_tpetra_block<myvec_t>::value, "");
  using sc_t = typename fix_t::ST;

  myvec_t B( *x_ );
  B.putScalar(1.22);

  myvec_t a(B);
  // a is wrapper of a block vector, so get a tpetra vector
  auto a_tp = a.data()->getVectorView();
  // now that we have a regular tpetra vector, we can check data
  Teuchos::ArrayRCP<const sc_t> dd = a_tp.getData();
  for (int i=0; i<a.localSize(); i++){
    EXPECT_DOUBLE_EQ( dd[i], 1.22 );
  }
}


TEST_F(tpetraBlockVectorGlobSize15BlockSize5Fixture,
       localSize){
  using namespace rompp;
  myvec_t v1( *x_ );
  EXPECT_EQ( numProc_, 3);
  EXPECT_EQ( v1.localSize(), 5);
}

TEST_F(tpetraBlockVectorGlobSize15BlockSize5Fixture,
       globalSize){
  using namespace rompp;
  myvec_t v1( *x_ );
  EXPECT_EQ( numProc_, 3);
  EXPECT_EQ( v1.globalSize(), 15);
}

TEST_F(tpetraBlockVectorGlobSize15BlockSize5Fixture,
       blockSize){
  using namespace rompp;
  myvec_t v1( *x_ );
  EXPECT_EQ( v1.data()->getBlockSize(), 4);
}

TEST_F(tpetraBlockVectorGlobSize15BlockSize5Fixture,
       SetScalar){
  using namespace rompp;
  using sc_t = typename fix_t::ST;

  myvec_t v1( *x_ );
  v1.putScalar(43.3);

  // v1 is a wrapper of a block vector, so get a tpetra vector
  auto v1_tp = v1.data()->getVectorView();

  // now that we have a regular tpetra vector, we can check data
  Teuchos::ArrayRCP<const sc_t> dd = v1_tp.getData();
  for (int i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( dd[i], 43.3 );
  }
}

TEST_F(tpetraBlockVectorGlobSize15BlockSize5Fixture,
       SetZero){
  using namespace rompp;
  using sc_t = typename fix_t::ST;

  using myvec_t = containers::Vector<native_t>;
  myvec_t v1( *x_ );
  v1.setZero();

  // v1 is a wrapper of a block vector, so get a tpetra vector
  auto v1_tp = v1.data()->getVectorView();

  // now that we have a regular tpetra vector, we can check data
  Teuchos::ArrayRCP<const sc_t> dd = v1_tp.getData();
  for (int i=0; i<v1.localSize(); i++){
    EXPECT_DOUBLE_EQ( dd[i], 0.0 );
  }
}
