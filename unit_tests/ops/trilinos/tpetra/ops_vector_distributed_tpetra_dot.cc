
#include "tpetra_only_fixtures.hpp"
#include "pressio_ops.hpp"

TEST_F(tpetraVectorGlobSize15Fixture, dot){
  using namespace pressio;

  using wrap_vec = typename tpetraVectorGlobSize15Fixture::vec_t;
  using myvec_t = containers::Vector<wrap_vec>;
  myvec_t a( *x_ );
  myvec_t b( *x_ );

  pressio::ops::fill(a, 1.0);
  pressio::ops::fill(b, 1.0);

  auto res = ::pressio::ops::dot(a, b);
  EXPECT_DOUBLE_EQ(res, 15.0);

  res = 0.0;
  ::pressio::ops::dot(a, b, res);
  EXPECT_DOUBLE_EQ(res, 15.0);
}
