
#include "../qr_utest_fixtures.hpp"
#include "QR_BASIC"

TEST_F(eigenDenseR9Fixture,
     HouseholderEigenDenseOutOfPlaceAndSolveEigenVecDynamic){
  using namespace rompp;

  fillMatrix();
  fillVector();

  //-------------------------------------------
  // R_type == void, in_place = false
  //-------------------------------------------
  using qr_algo = qr::Householder;
  qr::QRSolver<mymat_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  //  do Q^T * v_, i.e. project v_ onto Q
  myvec_t rhs(rompp::qr::test::numVectors_);
  qrObj.project(*v_, rhs);
  std::cout << " RHS" << std::setprecision(14)
  	    << *rhs.data() << std::endl;

  myvec_t y(rompp::qr::test::numVectors_);
  qrObj.solve(rhs, y);
  std::cout << std::setprecision(14)
  	    << " Y: "
  	    << *y.data() << std::endl;

  gold_.checkYForRsolve(y);
}
