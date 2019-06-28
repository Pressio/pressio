
#include "../qr_utest_fixtures.hpp"
#include "QR_BASIC"

TEST_F(epetraR9Fixture,
       TSQREpetraMVOutOfPlaceAndSolveEigenVecDynamic){
  using namespace rompp;

  using nat_v = Eigen::VectorXd;
  using myv_t = algebra::Vector<nat_v>;

  // default: R_type == void, in_place = false
  using qr_algo = qr::TSQR;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  //  do Q^T * v_, i.e. project v_ onto Q
  myv_t rhs(rompp::qr::test::numVectors_);
  qrObj.project(*v_, rhs);
  if (rank_==0)
    std::cout << " RHS" << std::setprecision(14)
	      << *rhs.data() << std::endl;

  // solve
  myv_t y(rompp::qr::test::numVectors_);
  qrObj.solve(rhs, y);
  if (rank_==0)
    std::cout << std::setprecision(14)
  	      << " Y: "
  	      << *y.data() << std::endl;

  gold_.checkYForRsolve(y);
}


TEST_F(epetraR9Fixture,
       TSQREpetraMVOutOfPlaceAndSolveEigenVecStatic){
  using namespace rompp;

  using nat_v = Eigen::Matrix<double,
			      rompp::qr::test::numVectors_, 1>;
  using myv_t = algebra::Vector<nat_v>;

  // default: R_type == void, in_place = false
  using qr_algo = qr::TSQR;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  //  do Q^T * v_, i.e. project v_ onto Q
  myv_t rhs;
  qrObj.project(*v_, rhs);
  if (rank_==0)
    std::cout << " RHS" << std::setprecision(14)
	      << *rhs.data() << std::endl;

  // solve
  myv_t y;
  qrObj.solve(rhs, y);
  if (rank_==0)
    std::cout << std::setprecision(14)
  	      << " Y: "
  	      << *y.data() << std::endl;

  gold_.checkYForRsolve(y);
}


TEST_F(epetraR9Fixture,
       TSQREpetraMVOutOfPlaceAndSolveTeuchosVector){
  using namespace rompp;

  using nat_v = Teuchos::SerialDenseVector<int,double>;
  using myv_t = algebra::Vector<nat_v>;

  //default: R_type == void, in_place = false
  using qr_algo = qr::TSQR;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  //  do Q^T * v_, i.e. project v_ onto Q
  myv_t rhs(rompp::qr::test::numVectors_);
  qrObj.project(*v_, rhs);
  if (rank_==0)
    std::cout << " RHS" << std::setprecision(14)
	      << *rhs.data() << std::endl;

  // solve
  myv_t y;
  y.resize(rompp::qr::test::numVectors_);
  qrObj.solve(rhs, y);
  if (rank_==0)
    std::cout << std::setprecision(14)
  	      << " Y: "
  	      << *y.data() << std::endl;

  gold_.checkYForRsolve(y);
}
