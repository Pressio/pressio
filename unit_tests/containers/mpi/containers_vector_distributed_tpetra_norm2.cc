
#include "tpetra_only_fixtures.hpp"

TEST_F(tpetraVectorGlobSize15Fixture, norm2){
  using namespace pressio;

  using wrap_vec = typename tpetraVectorGlobSize15Fixture::vec_t;
  using myvec_t = containers::Vector<wrap_vec>;
  myvec_t a( *x_ );
  a.putScalar(1.0);
  auto mynorm = containers::ops::norm2(a);
  EXPECT_DOUBLE_EQ(mynorm, std::sqrt(15.0));
}
