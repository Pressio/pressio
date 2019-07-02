
#include "block_tpetra_only_fixtures.hpp"

using fix_t	= tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture;
using native_t  = typename fix_t::mvec_t;
using mymvec_t	= pressio::containers::MultiVector<native_t>;

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       QueryWrappedData){
  using namespace pressio;

  mymvec_t A( *mv_ );
  ::testing::StaticAssertTypeEq<decltype(A.data()),
  				native_t * >();
  const mymvec_t A2( *mv_ );
  ::testing::StaticAssertTypeEq< decltype(A2.data()),
  				 const native_t * >();
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       Constructor){
  using namespace pressio;
  mymvec_t a( *mv_ );
}


TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       ConstructorFromNative){
  using namespace pressio;

  using sc_t = typename fix_t::ST;

  // set all values to 1.22
  mv_->putScalar(1.22);
  // create a wrapper
  mymvec_t A( *mv_ );

  // check that all values of A are 1.22
  // A is wrapper of a block MV, so get a tpetra MV
  auto Amv = A.data()->getMultiVectorView();

  // now that we have a regular tpetra MV, we can check data
  for (int j=0; j<A.localNumVectors(); j++){
    // get data for j column
    auto jCol_d = Amv.getData(j);
    // loop over rows and check
    for (int i=0; i<A.localLength(); i++){
      EXPECT_DOUBLE_EQ( jCol_d[i], 1.22 );
    }
  }
}


TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       CopyConstructor){
  using namespace pressio;

  using sc_t = typename fix_t::ST;

  // create a wrapper
  mymvec_t A( *mv_ );
  // set all values to 1.22
  A.data()->putScalar(1.22);

  // copy into B
  mymvec_t B( A );

  // check that all values of B are 1.22
  // B is wrapper of a block MV, so get a tpetra MV
  auto Bmv = B.data()->getMultiVectorView();

  // now that we have a regular tpetra MV, we can check data
  for (int j=0; j<B.localNumVectors(); j++){
    // get data for j column
    auto jCol_d = Bmv.getData(j);
    // loop over rows and check
    for (int i=0; i<B.localLength(); i++){
      EXPECT_DOUBLE_EQ( jCol_d[i], 1.22 );
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       localSize){
  using namespace pressio;
  mymvec_t A( *mv_ );

  // 5 is local length only if this is run with 3 proc
  EXPECT_EQ( numProc_, 3);
  EXPECT_EQ( A.localLength(), 5);
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       globalSize){
  using namespace pressio;
  mymvec_t A( *mv_ );
  // 5 is local length only if this is run with 3 proc
  EXPECT_EQ( numProc_, 3);
  EXPECT_EQ( A.globalLength(), 15);
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       localNumVecs){
  using namespace pressio;
  mymvec_t v1( *mv_ );
  EXPECT_EQ(v1.localNumVectors(), 3);
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       globalNumVecs){
  using namespace pressio;
  mymvec_t v1( *mv_ );
  EXPECT_EQ(v1.globalNumVectors(), 3);
}


TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       getMap){
  using namespace pressio;
  mymvec_t v1( *mv_ );
  auto const & mapO = v1.getDataMap();
  EXPECT_TRUE( mapO.isSameAs(*mv_->getMap()) );

  ::testing::StaticAssertTypeEq
      <decltype(mapO),
       const typename fix_t::map_t & >();
  EXPECT_TRUE(mapO.isContiguous());
}


TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       setZero){
  using namespace pressio;

  // create a wrapper
  mymvec_t B( *mv_ );
  // set all values to 1.22
  B.data()->putScalar(1.22);
  B.setZero();

  // check that all values are 0
  // B is wrapper of a block MV, so get a tpetra MV
  auto Bmv = B.data()->getMultiVectorView();

  // now that we have a regular tpetra MV, we can check data
  for (int j=0; j<B.localNumVectors(); j++){
    // get data for j column
    auto jCol_d = Bmv.getData(j);
    // loop over rows and check
    for (int i=0; i<B.localLength(); i++){
      EXPECT_DOUBLE_EQ( jCol_d[i], 0.0 );
    }
  }
}
