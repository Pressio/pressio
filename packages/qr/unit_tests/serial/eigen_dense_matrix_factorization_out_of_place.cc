
#include "../qr_utest_fixtures.hpp"
#include "QR_BASIC"

TEST_F(eigenDenseR9Fixture,
     HouseholderEigenDenseOutOfPlace){
  using namespace rompp;

  fillMatrix();

  //-------------------------------------------
  // R_type == void, in_place = false
  //-------------------------------------------
  using qr_algo = qr::Householder;
  qr::QRSolver<mymat_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  // Q is a core::MultiVector<...> by default
  const auto & Q = qrObj.cRefQFactor();
  EXPECT_EQ( Q.length(), qr::test::numRows_);
  EXPECT_EQ( Q.numVectors(), qr::test::numVectors_);
  std::cout << *Q.data() << "\n";
  // std::cout << *R.data() << "\n";

  for (auto i=0; i<qr::test::numRows_; i++){
    for (auto j=0; j<qr::test::numVectors_; j++){
      EXPECT_NEAR( std::abs(Q(i,j)),
  		   std::abs(gold_.trueQ_(i,j)), 1e-6);
    }
  }
}


TEST_F(eigenDenseR9Fixture,
       HouseholderEigenDenseOutOfPlaceStaticQMatrix){
  using namespace rompp;

  fillMatrix();

  //-------------------------------------------
  // R_type == void, in_place = false
  //-------------------------------------------
  using qr_algo = qr::Householder;
  qr::QRSolver<mymat_t, qr_algo, false, 4, 9, core::Matrix> qrObj;
  qrObj.computeThin( *A_ );

  // Q is a core::MultiVector<...> by default
  const auto & Q = qrObj.cRefQFactor();
  EXPECT_EQ( Q.rows(), qr::test::numRows_);
  EXPECT_EQ( Q.cols(), qr::test::numVectors_);
  std::cout << *Q.data() << "\n";

  for (auto i=0; i<qr::test::numRows_; i++){
    for (auto j=0; j<qr::test::numVectors_; j++){
      EXPECT_NEAR( std::abs(Q(i,j)),
  		   std::abs(gold_.trueQ_(i,j)), 1e-6);
    }
  }
}
