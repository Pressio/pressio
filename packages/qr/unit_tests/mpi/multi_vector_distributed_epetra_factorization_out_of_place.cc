
#include "../qr_utest_fixtures.hpp"
#include "QR_BASIC"

TEST_F(epetraR9Fixture,
       HouseholderEpetraMultiVectorOutOfPlace){
  using namespace rompp;

  // default: R_type == void, in_place = false
  using qr_algo = qr::Householder;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );
  const auto & Q = qrObj.cRefQFactor();
  Q.data()->Print(std::cout);
  checkQFactor(Q);
}

#ifdef HAVE_ANASAZI_TSQR
TEST_F(epetraR9Fixture,
       TSQRepetraMultiVectorOutOfPlace){
  using namespace rompp;

  // default: R_type == void, in_place = false
  using qr_algo = qr::TSQR;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );
  const auto & Q = qrObj.cRefQFactor();
  checkQFactor(Q);
}
#endif

#ifdef HAVE_BELOS_TSQR
TEST_F(epetraR9Fixture,
       BelosTSQRepetraMultiVectorOutOfPlace){
  using namespace rompp;

  // default: R_type == void, in_place = false
  using qr_algo = qr::TSQRBelos;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );
  const auto & Q = qrObj.cRefQFactor();
  checkQFactor(Q);
}
#endif

TEST_F(epetraR9Fixture,
       ModGrShEpetraMultiVectorOutOfPlace){
  using namespace rompp;

  // default: R_type == void, in_place = false
  using qr_algo = qr::ModifiedGramSchmidt;
  qr::QRSolver<mymvec_t, qr_algo> qrObj;
  qrObj.computeThin( *A_ );
  const auto & Q = qrObj.cRefQFactor();
  Q.data()->Print(std::cout);
  checkQFactor(Q);
}









// TEST_F(epetraR9Fixture,
//        TSQRepetraMultiVectorOutOfPlaceWrapREigen){
//   using namespace rompp;

//   // fill fixture matrix
//   fillMatrix();

//   //-------------------------------------------
//   // R_type == eigen_wrapper, in_place = false
//   //-------------------------------------------
//   using qr_algo = qr::TSQR;
//   using eig_mat = Eigen::Matrix<double,
// 				qr::test::numVectors_,
// 				qr::test::numVectors_>;
//   using R_type = core::Matrix<eig_mat>;
//   qr::QRSolverWrapR<mymvec_t, qr_algo, R_type> qrObj;
//   qrObj.computeThin( *A_ );

//   const auto & Q = qrObj.cRefQFactor();
//   EXPECT_EQ( Q.globalLength(), qr::test::numRows_);
//   EXPECT_EQ( Q.globalNumVectors(), qr::test::numVectors_);
//   for (auto i=0; i<localSize_; i++)
//     for (auto j=0; j<Q.localNumVectors(); j++){
//       EXPECT_NEAR( std::abs(Q(i,j)),
//   		   std::abs(gold_.trueQ_(i+shift_,j)), 1e-6);
//   }

//   // check R factor
//   const auto & R = qrObj.cRefRFactor();
//   if (rank_==0)  std::cout << "\n" << *R.data();
//   gold_.checkR(R);
// }
