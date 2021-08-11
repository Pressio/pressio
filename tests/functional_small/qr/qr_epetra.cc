
#include "fixtures.hpp"
#include "pressio/qr.hpp"

TEST_F(epetraR9Fixture, HouseholderEpetraMultiVectorOutOfPlace)
{
  using namespace pressio;

  // default: R_type == void, in_place = false
  using qr_algo = qr::Householder;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );
  const auto & Q = qrObj.cRefQFactor();
  Q.Print(std::cout);
  checkQFactor(Q);
}

TEST_F(epetraR9Fixture, TSQRepetraMultiVectorOutOfPlace)
{
  using namespace pressio;

  // default: R_type == void, in_place = false
  using qr_algo = qr::TSQR;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );
  const auto & Q = qrObj.cRefQFactor();
  Q.Print(std::cout);
  checkQFactor(Q);
}

TEST_F(epetraR9Fixture,
       ModGrShEpetraMultiVectorOutOfPlace)
{
  using namespace pressio;

  // default: R_type == void, in_place = false
  using qr_algo = qr::ModifiedGramSchmidt;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );
  const auto & Q = qrObj.cRefQFactor();
  Q.Print(std::cout);
  checkQFactor(Q);
}

TEST_F(epetraR9Fixture,
       TSQREpetraMVOutOfPlaceAndSolveTeuchosVector)
{
  using namespace pressio;

  //default: R_type == void, in_place = false
  using qr_algo = qr::TSQR;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  //  do Q^T * v_
  Teuchos::SerialDenseVector<int,double> rhs(pressio::qr::test::numVectors_);
  qrObj.applyQTranspose(*v_, rhs);

  // solve
  Teuchos::SerialDenseVector<int,double> y;
  y.resize(pressio::qr::test::numVectors_);
  qrObj.solve(rhs, y);

  EXPECT_EQ(y.length(), 4);
  gold_.checkYForRsolve(y);
}

#ifdef PRESSIO_ENABLE_TPL_EIGEN
TEST_F(epetraR9Fixture,
       TSQREpetraMVOutOfPlaceAndSolveEigenVecDynamic)
{
  using namespace pressio;

  // default: R_type == void, in_place = false
  using qr_algo = qr::TSQR;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin(*A_);

  //  do Q^T * v_
  Eigen::VectorXd rhs(pressio::qr::test::numVectors_);
  qrObj.applyQTranspose(*v_, rhs);
  if (rank_==0)
    std::cout << " RHS" << std::setprecision(14) << rhs << std::endl;

  // solve
  Eigen::VectorXd y(pressio::qr::test::numVectors_);
  qrObj.solve(rhs, y);
  if (rank_==0)
    std::cout << std::setprecision(14) << " Y: " << y << std::endl;

  EXPECT_EQ(y.size(), 4);
  gold_.checkYForRsolve(y);
}

TEST_F(epetraR9Fixture,
       TSQREpetraMVOutOfPlaceAndSolveEigenVecStatic)
{
  using namespace pressio;

  using myv_t = Eigen::Matrix<double, pressio::qr::test::numVectors_, 1>;

  // default: R_type == void, in_place = false
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
