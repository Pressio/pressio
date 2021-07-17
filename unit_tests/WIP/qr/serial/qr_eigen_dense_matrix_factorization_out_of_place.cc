
#include "../qr_utest_fixtures.hpp"
#include "pressio_qr.hpp"

TEST_F(eigenDenseR9Fixture, HouseholderEigenDenseOutOfPlace){
  using namespace pressio;

  // defaults to: R_type == void, in_place = false
  using qr_algo = qr::Householder;
  qr::QRSolver<matrix_type, qr_algo> qrObj;
  qrObj.computeThin( A_ );

  const auto & Q = qrObj.cRefQFactor();
  checkQFactor(Q);
}
