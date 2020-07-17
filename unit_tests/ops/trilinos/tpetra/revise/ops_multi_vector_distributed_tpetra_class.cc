
#include "tpetra_only_fixtures.hpp"

TEST_F(tpetraMultiVectorGlobSize15Fixture,
       SetZero){
  using namespace pressio;
  using sc_t = typename tpetraMultiVectorGlobSize15Fixture::ST;

  using mymvec_t = containers::MultiVector<typename tpetraMultiVectorGlobSize15Fixture::mvec_t>;
  mymvec_t v1( *x_ );
  ::pressio::ops::set_zero(v1);

  for (int k=0; k<4; k++){
    Teuchos::ArrayRCP<const sc_t> dd = v1.data()->getData(k);
    for (int i=0; i<v1.extentLocal(0); i++){
      EXPECT_DOUBLE_EQ( dd[i], 0.0 );
    }
  }
}