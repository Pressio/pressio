
#include "epetra_only_fixtures.hpp"
#include "pressio_ops.hpp"

TEST_F(epetraVectorGlobSize15Fixture, norm2){
  using namespace pressio;

  using myvec_t = containers::Vector<Epetra_Vector>;
  myvec_t a( *contigMap_ );
  a.data()->PutScalar(1.0);
  auto mynorm = ops::norm2(a);
  EXPECT_DOUBLE_EQ(mynorm, std::sqrt(15.0));
}
