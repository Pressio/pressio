
#include "../qr_utest_fixtures.hpp"
#include "pressio_qr.hpp"

TEST_F(epetraR9Fixture,
       TSQREpetraMVOutOfPlaceAndSolveEigenVecDynamic){
  using namespace pressio;

  using nat_v = Eigen::VectorXd;
  using myv_t = containers::Vector<nat_v>;

  // default: R_type == void, in_place = false
  using qr_algo = qr::TSQR;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  //  do Q^T * v_, i.e. project v_ onto Q
  myv_t rhs(pressio::qr::test::numVectors_);
  qrObj.applyQTranspose(*v_, rhs);
  if (rank_==0)
    std::cout << " RHS" << std::setprecision(14)
	      << *rhs.data() << std::endl;

  // solve
  myv_t y(pressio::qr::test::numVectors_);
  qrObj.solve(rhs, y);
  if (rank_==0)
    std::cout << std::setprecision(14)
  	      << " Y: "
  	      << *y.data() << std::endl;

  gold_.checkYForRsolve(y);
}


TEST_F(epetraR9Fixture,
       TSQREpetraMVOutOfPlaceAndSolveEigenVecStatic){
  using namespace pressio;

  using nat_v = Eigen::Matrix<double, pressio::qr::test::numVectors_, 1>;
  using myv_t = containers::Vector<nat_v>;

  // default: R_type == void, in_place = false
  using qr_algo = qr::TSQR;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  //  do Q^T * v_, i.e. project v_ onto Q
  myv_t rhs;
  qrObj.applyQTranspose(*v_, rhs);
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
  using namespace pressio;

  using nat_v = Teuchos::SerialDenseVector<int,double>;
  using myv_t = containers::Vector<nat_v>;

  //default: R_type == void, in_place = false
  using qr_algo = qr::TSQR;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  //  do Q^T * v_, i.e. project v_ onto Q
  myv_t rhs(pressio::qr::test::numVectors_);
  qrObj.applyQTranspose(*v_, rhs);
  // if (rank_==0)
  //   std::cout << " RHS" << std::setprecision(14)
	 //      << *rhs.data() << std::endl;

  // solve
  myv_t y;
  y.data()->resize(pressio::qr::test::numVectors_);
  qrObj.solve(rhs, y);
  // if (rank_==0)
  //   std::cout << std::setprecision(14)
  // 	      << " Y: "
  // 	      << *y.data() << std::endl;

  gold_.checkYForRsolve(y);
}
