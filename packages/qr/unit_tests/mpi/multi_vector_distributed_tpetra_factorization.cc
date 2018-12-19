
#include "../qr_utest_fixtures.hpp"
#include "../../src/qr_tpetra_multi_vector_householder.hpp"
#include "../../src/qr_tpetra_multi_vector_tsqr.hpp"


TEST_F(tpetraR9Fixture,
       TpetraMultiVectorQRFactorizationHouseholder){
  using namespace rompp;

  // fill matrix inside fixture
  fillMatrix();

  // do QR
  using qr_algo = qr::Householder;
  using R_type = core::Matrix<Eigen::MatrixXd>;
  qr::QRSolver<mymvec_t, R_type, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  const auto & Q = qrObj.cRefQFactor();
  EXPECT_EQ( Q.globalLength(), 9);
  EXPECT_EQ( Q.globalNumVectors(), 4);

  int localSize = (rank_==0) ? 5 : 4;
  int shift = (rank_==0) ? 0 : 5;
  for (auto j=0; j<Q.localNumVectors(); j++){
    auto colData = Q.data()->getData(j);
    for (auto i=0; i<localSize; i++){
      EXPECT_NEAR( std::abs(colData[i]),
		   std::abs(gold.trueQ_(i+shift,j)), 1e-6);
    }
  }

  // check R factor
  const auto & R = qrObj.cRefRFactor();
  EXPECT_EQ( R.cols(), 4 );
  if (rank_==0)  std::cout << "\n" << *R.data();

  for (auto i=0; i<4; i++)
    for (auto j=0; j<4; j++)
      EXPECT_NEAR( std::abs(R(i,j)), std::abs(gold.trueR_(i,j)), 1e-6);
}



TEST_F(tpetraR9Fixture,
       TpetraMultiVectorQRFactorizationTSQR){
  using namespace rompp;

  // fill matrix inside fixture
  fillMatrix();

  // do QR
  using qr_algo = qr::TSQR;
  constexpr int nC = 4;
  assert(nC == numVectors_);
  using R_type = rompp::core::Matrix<Eigen::Matrix<double, nC, nC>>;

  qr::QRSolver<mymvec_t, R_type, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  const auto & Q = qrObj.cRefQFactor();
  EXPECT_EQ( Q.globalLength(), 9);
  EXPECT_EQ( Q.globalNumVectors(), 4);

  int localSize = (rank_==0) ? 5 : 4;
  int shift = (rank_==0) ? 0 : 5;
  for (auto j=0; j<Q.localNumVectors(); j++){
    auto colData = Q.data()->getData(j);
    for (auto i=0; i<localSize; i++){
      EXPECT_NEAR( std::abs(colData[i]),
		   std::abs(gold.trueQ_(i+shift,j)), 1e-6);
    }
  }

  // check R factor
  const auto & R = qrObj.cRefRFactor();
  EXPECT_EQ( R.cols(), 4 );
  if (rank_==0)  std::cout << "\n" << *R.data();

  for (auto i=0; i<4; i++)
    for (auto j=0; j<4; j++)
      EXPECT_NEAR( std::abs(R(i,j)), std::abs(gold.trueR_(i,j)), 1e-6);
}
