
#include "epetra_only_fixtures.hpp"

TEST_F(epetraVectorGlobSize15Fixture, norm2){
  using namespace rompp;

  using myvec_t = algebra::Vector<Epetra_Vector>;
  myvec_t a( *contigMap_ );
  a.putScalar(1.0);
  auto mynorm = algebra::ops::norm2(a);
  EXPECT_DOUBLE_EQ(mynorm, std::sqrt(15.0));
}
