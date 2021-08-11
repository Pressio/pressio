
#include "fixtures.hpp"
#include "pressio/qr.hpp"

TEST_F(tpetraR9Fixture,
       HouseholderTpetraMultiVectorOutOfPlace)
{
  using namespace pressio;

  // default: R_type == void, in_place = false
  using qr_algo = qr::Householder;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  const auto & Q = qrObj.cRefQFactor();
  checkQFactor(Q);
}

TEST_F(tpetraR9Fixture,
       TSQRtpetraMultiVectorOutOfPlace)
{
  using namespace pressio;

  // default: R_type == void, in_place = false
  using qr_algo = qr::TSQR;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  const auto & Q = qrObj.cRefQFactor();
  checkQFactor(Q);
}

TEST_F(tpetraR9Fixture,
       ModGrShTpetraMultiVectorOutOfPlace)
{
  using namespace pressio;

  // default: R_type == void, in_place = false
  using qr_algo = qr::ModifiedGramSchmidt;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  const auto & Q = qrObj.cRefQFactor();
  checkQFactor(Q);
}

#ifdef PRESSIO_ENABLE_TPL_EIGEN
TEST_F(tpetraR9Fixture,
       TSQRTpetraMVOutOfPlaceAndSolveEigenVecDynamic)
{
  using namespace pressio;

  // default: R_type == void, in_place = false
  using qr_algo = qr::TSQR;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  //  do Q^T * v_
  Eigen::VectorXd rhs(pressio::qr::test::numVectors_);
  qrObj.applyQTranspose(*v_, rhs);
  if (rank_==0)
    std::cout << " RHS" << std::setprecision(14)
        << *rhs.data() << std::endl;

  // solve
  Eigen::VectorXd y(pressio::qr::test::numVectors_);
  qrObj.solve(rhs, y);
  if (rank_==0)
    std::cout << std::setprecision(14)
          << " Y: "
          << *y.data() << std::endl;

  gold_.checkYForRsolve(y);
}

TEST_F(tpetraR9Fixture,
       TSQRTpetraMVOutOfPlaceAndSolveEigenVecStatic)
{
  using namespace pressio;
  using myv_t = Eigen::Matrix<double, qr::test::numVectors_, 1>;

  // detault: R_type == void, in_place = false
  using qr_algo = qr::TSQR;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  //  do Q^T * v_
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
#endif

