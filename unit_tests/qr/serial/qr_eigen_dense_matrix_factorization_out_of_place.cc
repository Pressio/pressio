
#include "../qr_utest_fixtures.hpp"
#include "pressio_qr.hpp"

TEST_F(eigenDenseR9Fixture, HouseholderEigenDenseOutOfPlace){
  using namespace pressio;

  // R_type == void, in_place = false
  using qr_algo = qr::Householder;
  qr::QRSolver<mymat_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  // Q is a containers::MultiVector<...> by default
  const auto & Q = qrObj.cRefQFactor();
  checkQFactor(*Q.data());
}
