
#include "../qr_utest_fixtures.hpp"
#include "../../src/qr_tpetra_multi_vector_householder.hpp"
#include "../../src/qr_tpetra_multi_vector_tsqr.hpp"


TEST_F(tpetraR9Fixture,
       TpetraMultiVectorTSQRprojectDynamicEigenVecResult){
  using namespace rompp;

  using eig_col_vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  // fill matrix inside fixture
  fillMatrix();
  fillVector();

  constexpr int nC = 4;
  assert(nC == numVectors_);
  using qr_algo = qr::TSQR;
  using R_type = containers::Matrix<Eigen::Matrix<double, nC, nC>>;
  qr::QRSolver<mymvec_t, R_type, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  //  do Q^T * v_ = y, i.e. project v_ onto Q
  containers::Vector<eig_col_vec> rhs;
  qrObj.project(*v_, rhs);
  if (rank_==0)
    std::cout << std::setprecision(14) << *rhs.data() << std::endl;

  // gold solution
  rompp::qr::test::qrGoldSol<double> gold;

  // solve
  containers::Vector<eig_col_vec> y;
  qrObj.solve(rhs, y);
  if (rank_==0)
    std::cout << std::setprecision(14) << *y.data() << std::endl;

  EXPECT_EQ(y.size(), 4);
  for (auto i=0; i<y.size(); i++)
    EXPECT_NEAR( std::abs(gold.trueYForRSolve_[i]),
  		 std::abs(y[i]), 1e-12 );
}


TEST_F(tpetraR9Fixture,
       TpetraMultiVectorTSQRprojectStaticEigenVecResult){
  using namespace rompp;

  using eig_col_vec = Eigen::Matrix<double, 4, 1>;

  // fill matrix inside fixture
  fillMatrix();
  fillVector();

  constexpr int nC = 4;
  assert(nC == numVectors_);
  using qr_algo = qr::TSQR;
  using R_type = containers::Matrix<Eigen::Matrix<double, nC, nC>>;
  qr::QRSolver<mymvec_t, R_type, qr_algo> qrObj;
  qrObj.computeThin( *A_ );

  //  do Q^T * v_ = y, i.e. project v_ onto Q
  containers::Vector<eig_col_vec> rhs;
  qrObj.project(*v_, rhs);
  if (rank_==0)
    std::cout << std::setprecision(14) << *rhs.data() << std::endl;

  // gold solution
  rompp::qr::test::qrGoldSol<double> gold;

  // solve
  containers::Vector<eig_col_vec> y;
  qrObj.solve(rhs, y);
  if (rank_==0)
    std::cout << std::setprecision(14) << *y.data() << std::endl;

  EXPECT_EQ(y.size(), 4);
  for (auto i=0; i<y.size(); i++)
    EXPECT_NEAR( std::abs(gold.trueYForRSolve_[i]),
  		 std::abs(y[i]), 1e-12 );
}
