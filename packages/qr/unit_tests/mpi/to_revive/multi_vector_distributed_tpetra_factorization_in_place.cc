
#include "../qr_utest_fixtures.hpp"
#include "QR_BASIC"

TEST_F(tpetraR9Fixture,
       TSQRtpetraMultiVectorInPlace){
  using namespace pressio;

  using qr_algo = qr::TSQR;
  constexpr bool in_place = true;

  //-------------------------------------------
  // R_type == void, in_place = true
  //-------------------------------------------
  qr::QRSolver<mymvec_t, qr_algo, in_place> qrObj;
  qrObj.computeThin( *A_ );

  EXPECT_EQ( A_->globalLength(), qr::test::numRows_);
  EXPECT_EQ( A_->globalNumVectors(), qr::test::numVectors_);

  for (auto j=0; j<A_->localNumVectors(); j++){
    auto colData = A_->data()->getData(j);
    for (auto i=0; i<localSize_; i++){
      EXPECT_NEAR( std::abs(colData[i]),
                   std::abs(gold_.trueQ_(i+shift_,j)),
                   1e-6);
    }
  }
}



TEST_F(tpetraR9Fixture,
       TSQRepetraMultiVectorInPlaceWrapREigen){
  using namespace pressio;

  using qr_algo = qr::TSQR;
  constexpr bool in_place = true;

  //-------------------------------------------
  // R_type == eigen_wrapper, in_place = true
  //-------------------------------------------
  using eig_mat = Eigen::Matrix<double,
				qr::test::numVectors_,
				qr::test::numVectors_>;
  using R_type = containers::Matrix<eig_mat>;
  qr::QRSolverWrapR<mymvec_t, qr_algo, R_type, in_place> qrObj;
  qrObj.computeThin( *A_ );

  EXPECT_EQ( A_->globalLength(), qr::test::numRows_);
  EXPECT_EQ( A_->globalNumVectors(), qr::test::numVectors_);
  for (auto j=0; j<A_->localNumVectors(); j++){
    auto colData = A_->data()->getData(j);
    for (auto i=0; i<localSize_; i++){
      EXPECT_NEAR( std::abs(colData[i]),
                   std::abs(gold_.trueQ_(i+shift_,j)),
                   1e-6);
    }
  }

  // check R factor
  const auto & R = qrObj.cRefRFactor();
  if (rank_==0)  std::cout << "\n" << *R.data();
  gold_.checkR(R);
}



TEST_F(tpetraR9Fixture,
       TSQRepetraMultiVectorInPlaceWrapRTeuchosMatrix){
  using namespace pressio;

  using qr_algo = qr::TSQR;
  constexpr bool in_place = true;

  //-------------------------------------------
  // R_type == eigen_wrapper, in_place = true
  //-------------------------------------------
  using nat_r_mat = Teuchos::SerialDenseMatrix<int, double>;
  using R_type = containers::Matrix<nat_r_mat>;
  qr::QRSolverWrapR<mymvec_t, qr_algo, R_type, in_place> qrObj;
  qrObj.computeThin( *A_ );

  EXPECT_EQ( A_->globalLength(), qr::test::numRows_);
  EXPECT_EQ( A_->globalNumVectors(), qr::test::numVectors_);
  for (auto j=0; j<A_->localNumVectors(); j++){
    auto colData = A_->data()->getData(j);
    for (auto i=0; i<localSize_; i++){
      EXPECT_NEAR( std::abs(colData[i]),
                   std::abs(gold_.trueQ_(i+shift_,j)),
                   1e-6);
    }
  }

  // check R factor
  const auto & R = qrObj.cRefRFactor();
  if (rank_==0)  std::cout << "\n" << *R.data();
  gold_.checkR(R);
}















// TEST_F(epetraR9Fixture,
//        QRepetraMultiVectorTSQRStoreTeuchosRmat){
//   using namespace pressio;

//   // fill matrix inside fixture
//   fillMatrix();

//   constexpr int nC = 4;
//   assert(nC == numVectors_);

//   // do QR
//   using qr_algo = qr::TSQR;
//   using nat_r_mat = Teuchos::SerialDenseMatrix<int, double>;
//   using R_type = containers::Matrix<nat_r_mat>;
//   qr::QRSolver<mymvec_t, R_type, qr_algo> qrObj;

//   //  compute
//   qrObj.computeThin( *A_ );

//   const auto & Q = qrObj.cRefQFactor();
//   EXPECT_EQ( Q.globalLength(), 9);
//   EXPECT_EQ( Q.globalNumVectors(), 4);

//   // gold solution
//   pressio::qr::test::qrGoldSol<double> gold;

//   int localSize = (rank_==0) ? 5 : 4;
//   int shift = (rank_==0) ? 0 : 5;
//   for (auto i=0; i<localSize; i++)
//     for (auto j=0; j<Q.localNumVectors(); j++){
//       EXPECT_NEAR( std::abs(Q(i,j)), std::abs(gold.trueQ_(i+shift,j)), 1e-6);
//   }

//   // check R factor
//   const auto & R = qrObj.cRefRFactor();
//   EXPECT_EQ( R.cols(), 4 );
//   EXPECT_EQ( R.rows(), 4 );
//   if (rank_==0)  std::cout << "\n" << *R.data();

//   for (auto i=0; i<4; i++)
//     for (auto j=0; j<4; j++)
//       EXPECT_NEAR( std::abs(R(i,j)), std::abs(gold.trueR_(i,j)), 1e-6);
// }
