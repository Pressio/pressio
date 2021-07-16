
#include "../qr_utest_fixtures.hpp"
#include "pressio_qr.hpp"

TEST_F(eigenDenseR9Fixture,
     HouseholderEigenDenseOutOfPlaceAndSolveEigenVecDynamic)
{
  using namespace pressio;

  // R_type == void, in_place = false
  using qr_algo = qr::Householder;
  qr::QRSolver<matrix_type, qr_algo> qrObj;
  qrObj.computeThin( A_ );

  vector_type rhs(pressio::qr::test::numVectors_);
  qrObj.applyQTranspose(v_, rhs);
  std::cout << " RHS" << std::setprecision(14)
  	    << rhs << std::endl;

  vector_type y(pressio::qr::test::numVectors_);
  qrObj.solve(rhs, y);
  std::cout << std::setprecision(14)
  	    << " Y: "
  	    << y << std::endl;

  gold_.checkYForRsolve(y);
}
