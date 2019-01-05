
#include "../qr_utest_fixtures.hpp"
#include "QR_BASIC"

TEST_F(eigenDenseR9Fixture,
     HouseholderEigenDenseOutOfPlace){
  using namespace rompp;

  // R_type == void, in_place = false
  using qr_algo = qr::Householder;
  qr::QRSolver<mymat_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  // Q is a core::MultiVector<...> by default
  const auto & Q = qrObj.cRefQFactor();
  checkQFactor(*Q.data());
}


TEST_F(eigenDenseR9Fixture,
       HouseholderEigenDenseOutOfPlaceStaticQMatrix){
  using namespace rompp;

  // R_type == void, in_place = false
  using qr_algo = qr::Householder;
  qr::QRSolver<mymat_t, qr_algo, false, 4, 9, core::Matrix> qrObj;
  qrObj.computeThin( *A_ );

  // Q is a core::MultiVector<...> by default
  const auto & Q = qrObj.cRefQFactor();
  checkQFactor(*Q.data());
}
