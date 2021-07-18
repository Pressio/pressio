
#include "epetra_only_fixtures.hpp"
#include "pressio_ops.hpp"

TEST_F(epetraVectorGlobSize15Fixture, dot){
  using namespace pressio;

  using myvec_t = containers::Vector<Epetra_Vector>;
  myvec_t a( *contigMap_ );
  myvec_t b( *contigMap_ );

  a.data()->PutScalar(1.0);
  b.data()->PutScalar(1.0);

  auto res = ::pressio::ops::dot(a, b);
  EXPECT_DOUBLE_EQ(res, 15.0);

  res = 0.0;
  ::pressio::ops::dot(a, b, res);
  EXPECT_DOUBLE_EQ(res, 15.0);
}
