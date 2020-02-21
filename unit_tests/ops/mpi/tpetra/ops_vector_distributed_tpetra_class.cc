
#include "tpetra_only_fixtures.hpp"
#include "pressio_ops.hpp"

using native_t = typename tpetraVectorGlobSize15Fixture::vec_t;


TEST_F(tpetraVectorGlobSize15Fixture,
       SetScalar){
  using namespace pressio;
  using sc_t = typename tpetraVectorGlobSize15Fixture::ST;

  using myvec_t = containers::Vector<native_t>;
  myvec_t v1( *x_ );
  pressio::ops::fill(v1, 43.3);

  Teuchos::ArrayRCP<const sc_t> dd = v1.data()->getData();
  for (int i=0; i<v1.extentLocal(0); i++){
    EXPECT_DOUBLE_EQ( dd[i], 43.3 );
  }
}

TEST_F(tpetraVectorGlobSize15Fixture,
       SetZero){
  using namespace pressio;
  using sc_t = typename tpetraVectorGlobSize15Fixture::ST;

  using myvec_t = containers::Vector<native_t>;
  myvec_t v1( *x_ );
  ::pressio::ops::set_zero(v1);

  Teuchos::ArrayRCP<const sc_t> dd = v1.data()->getData();
  for (int i=0; i<v1.extentLocal(0); i++){
    EXPECT_DOUBLE_EQ( dd[i], 0.0 );
  }
}
