
#include "../qr_utest_fixtures.hpp"
#include "../../src/qr_epetra_multi_vector_householder.hpp"
#include "../../src/qr_epetra_multi_vector_tsqr.hpp"


TEST_F(epetraR9Fixture,
       EpetraMultiVectorQRFactorizationHouseholder){
  using namespace rompp;

  // fill matrix inside fixture
  fillMatrix();

  // do QR
  using qr_algo = qr::Householder;
  using R_type = core::Matrix<Eigen::MatrixXd>;
  qr::QRSolver<mymvec_t, R_type, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  const auto & Q = qrObj.cRefQFactor();
  //Q.data()->Print(std::cout << std::setprecision(4));
  EXPECT_EQ( Q.globalLength(), 9);
  EXPECT_EQ( Q.globalNumVectors(), 4);

  // gold solution
  rompp::qr::test::qrGoldSol<double> gold;

  int localSize = (rank_==0) ? 5 : 4;
  int shift = (rank_==0) ? 0 : 5;
  for (auto i=0; i<localSize; i++)
    for (auto j=0; j<Q.localNumVectors(); j++){
      EXPECT_NEAR( std::abs(Q(i,j)), std::abs(gold.trueQ_(i+shift,j)), 1e-6);
  }

  // check R factor
  const auto & R = qrObj.cRefRFactor();
  EXPECT_EQ( R.cols(), 4 );
  if (rank_==0)  std::cout << std::setprecision(14) << *R.data();
  std::cout << std::endl;

  for (auto i=0; i<4; i++)
    for (auto j=0; j<4; j++)
      EXPECT_NEAR( std::abs(R(i,j)), std::abs(gold.trueR_(i,j)), 1e-6);
}



TEST_F(epetraR9Fixture,
       EpetraMultiVectorQRFactorizationTSQR){
  using namespace rompp;

  // fill matrix inside fixture
  fillMatrix();

  constexpr int nC = 4;
  assert(nC == numVectors_);

  // do QR
  using qr_algo = qr::TSQR;
  using R_type = core::Matrix<Eigen::Matrix<double, nC, nC>>;
  qr::QRSolver<mymvec_t, R_type, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  const auto & Q = qrObj.cRefQFactor();
  //Q.data()->Print(std::cout << std::setprecision(4));
  EXPECT_EQ( Q.globalLength(), 9);
  EXPECT_EQ( Q.globalNumVectors(), 4);

  // gold solution
  rompp::qr::test::qrGoldSol<double> gold;

  int localSize = (rank_==0) ? 5 : 4;
  int shift = (rank_==0) ? 0 : 5;
  for (auto i=0; i<localSize; i++)
    for (auto j=0; j<Q.localNumVectors(); j++){
      EXPECT_NEAR( std::abs(Q(i,j)), std::abs(gold.trueQ_(i+shift,j)), 1e-6);
  }

  // check R factor
  const auto & R = qrObj.cRefRFactor();
  EXPECT_EQ( R.cols(), 4 );
  EXPECT_EQ( R.rows(), 4 );
  if (rank_==0)  std::cout << "\n" << *R.data();

  for (auto i=0; i<4; i++)
    for (auto j=0; j<4; j++)
      EXPECT_NEAR( std::abs(R(i,j)), std::abs(gold.trueR_(i,j)), 1e-6);
}
