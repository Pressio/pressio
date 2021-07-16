
#include "tpetra_only_fixtures.hpp"

using native_t = typename tpetraVectorGlobSize15Fixture::vec_t;

TEST_F(tpetraVectorGlobSize15Fixture,
       QueryWrappedData){
  using namespace pressio;

  using myvec_t = containers::Vector<native_t>;
  myvec_t v1( *x_ );
  ::testing::StaticAssertTypeEq<decltype(v1.data()),
          native_t * >();
  const myvec_t v2( *x_ );
  ::testing::StaticAssertTypeEq< decltype(v2.data()),
           const native_t * >();
}

TEST_F(tpetraVectorGlobSize15Fixture,
       localSize){
  using namespace pressio;
  using myvec_t = containers::Vector<native_t>;
  myvec_t v1( *x_ );
  EXPECT_EQ( v1.extentLocal(0), 5);
}

TEST_F(tpetraVectorGlobSize15Fixture,
       globalSize){
  using namespace pressio;
  using myvec_t = containers::Vector<native_t>;
  myvec_t v1( *x_ );
  EXPECT_EQ( v1.extent(0), 15);
}

