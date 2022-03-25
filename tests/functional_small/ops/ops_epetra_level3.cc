
#include "epetra_only_fixtures.hpp"
#include "pressio/ops.hpp"

#ifdef PRESSIO_ENABLE_TPL_EIGEN
TEST_F(epetraMultiVectorGlobSize15Fixture, mv_T_mv_storein_eigen_C)
{
    auto A = pressio::ops::clone(*myMv_);
    std::array<double, 4> ac{1.,2.,3.,4.};
    for (int i=0; i<A.NumVectors(); ++i) {
      A(i)->PutScalar(ac[i]);
    }

    Epetra_MultiVector B(*contigMap_, 3);
    std::array<double, 3> bc{1.2, 2.2, 3.2};
    for (int i=0; i<3; ++i) {
      B(i)->PutScalar(bc[i]);
    }

    Eigen::MatrixXd C(A.NumVectors(), B.NumVectors());
    C.setConstant(0.);

    // C = 1*C + 1.5 A^T B
    pressio::ops::product(pressio::transpose(), pressio::nontranspose(),
        1.5, A, B, 1.0, C);

    auto C2 = pressio::ops::product<Eigen::MatrixXd>(pressio::transpose(), pressio::nontranspose(),
        1.5, A, B);

    for (auto i=0; i<C.rows(); i++){
        for (auto j=0; j<C.cols(); j++){
            const auto gold = ac[i]*A.GlobalLength()*1.5*bc[j];
            EXPECT_NEAR( C(i,j), gold, 1e-12);
            EXPECT_NEAR( C2(i,j), gold, 1e-12);
        }
    }
}

TEST_F(epetraMultiVectorGlobSize15Fixture, mv_T_mv_storein_eigen_C_beta0)
{
    auto A = pressio::ops::clone(*myMv_);
    std::array<double, 4> ac{1.,2.,3.,4.};
    for (int i=0; i<A.NumVectors(); ++i) {
      A(i)->PutScalar(ac[i]);
    }

    Epetra_MultiVector B(*contigMap_, 3);
    std::array<double, 3> bc{1.2, 2.2, 3.2};
    for (int i=0; i<3; ++i) {
      B(i)->PutScalar(bc[i]);
    }

    Eigen::MatrixXd C(A.NumVectors(), B.NumVectors());
    C.setConstant(NAN);

    // C = 0*NAN + 1.5 A^T B
    pressio::ops::product(pressio::transpose(), pressio::nontranspose(),
        1.5, A, B, 0.0, C);

    auto C2 = pressio::ops::product<Eigen::MatrixXd>(pressio::transpose(), pressio::nontranspose(),
        1.5, A, B);

    for (auto i=0; i<C.rows(); i++){
        for (auto j=0; j<C.cols(); j++){
            const auto gold = ac[i]*A.GlobalLength()*1.5*bc[j];
            EXPECT_NEAR( C(i,j), gold, 1e-12);
            EXPECT_NEAR( C2(i,j), gold, 1e-12);
        }
    }
}

TEST_F(epetraMultiVectorGlobSize15Fixture, mv_T_self_storein_eigen_C)
{
    auto A = pressio::ops::clone(*myMv_);
    std::array<double, 4> ac{1.,2.,3.,4.};
    for (int i=0; i<A.NumVectors(); ++i) {
      A(i)->PutScalar(ac[i]);
    }

    Eigen::MatrixXd C(A.NumVectors(), A.NumVectors());
    C.setConstant(0.);

    // C = 1*C + 1.5 A^T A
    pressio::ops::product(pressio::transpose(), pressio::nontranspose(),
        1.5, A, 1.0, C);

    if(rank_==0){
        std::cout << C << std::endl;
    }

    for (auto i=0; i<C.rows(); i++){
        for (auto j=0; j<C.cols(); j++){
            const auto gold = ac[i]*A.GlobalLength()*1.5*ac[j];
            EXPECT_NEAR( C(i,j), gold, 1e-12);
        }
    }
}

TEST_F(epetraMultiVectorGlobSize15Fixture, mv_T_self_storein_eigen_C_beta0)
{
    auto A = pressio::ops::clone(*myMv_);
    std::array<double, 4> ac{1.,2.,3.,4.};
    for (int i=0; i<A.NumVectors(); ++i) {
      A(i)->PutScalar(ac[i]);
    }

    Eigen::MatrixXd C(A.NumVectors(), A.NumVectors());
    C.setConstant(NAN);

    // C = 0*NAN + 1.5 A^T A
    pressio::ops::product(pressio::transpose(), pressio::nontranspose(),
        1.5, A, 0.0, C);

    if(rank_==0){
        std::cout << C << std::endl;
    }

    for (auto i=0; i<C.rows(); i++){
        for (auto j=0; j<C.cols(); j++){
            const auto gold = ac[i]*A.GlobalLength()*1.5*ac[j];
            EXPECT_NEAR( C(i,j), gold, 1e-12);
        }
    }
}
#endif
