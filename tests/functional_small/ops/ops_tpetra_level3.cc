
#include "tpetra_only_fixtures.hpp"
#include "pressio/ops.hpp"

#ifdef PRESSIO_ENABLE_TPL_EIGEN
TEST_F(tpetraMultiVectorGlobSize15Fixture, mv_T_mv_storein_eigen_C)
{
    auto A = pressio::ops::clone(*myMv_);
    std::array<double, 4> ac{1.,2.,3.,4.};
    for (std::size_t i=0; i<A.getNumVectors(); ++i) {
      A.getVectorNonConst(i)->putScalar(ac[i]);
    }

    mvec_t B(contigMap_, 3);
    std::array<double, 3> bc{1.2, 2.2, 3.2};
    for (int i=0; i<3; ++i) {
      B.getVectorNonConst(i)->putScalar(bc[i]);
    }

    Eigen::MatrixXd C(A.getNumVectors(), B.getNumVectors());
    C.setConstant(0.);

    // C = 1*C + 1.5 A^T B
    pressio::ops::product(
        pressio::transpose(),
        pressio::nontranspose(),
        1.5, A, B, 1.0, C);

    for (auto i=0; i<C.rows(); i++){
        for (auto j=0; j<C.cols(); j++){
            const auto gold = ac[i]*A.getGlobalLength()*1.5*bc[j];
            EXPECT_NEAR( C(i,j), gold, 1e-12);
        }
    }
}

TEST_F(tpetraMultiVectorGlobSize15Fixture, mv_T_mv_storein_eigen_C_beta0)
{
    auto A = pressio::ops::clone(*myMv_);
    std::array<double, 4> ac{1.,2.,3.,4.};
    for (std::size_t i=0; i<A.getNumVectors(); ++i) {
      A.getVectorNonConst(i)->putScalar(ac[i]);
    }

    mvec_t B(contigMap_, 3);
    std::array<double, 3> bc{1.2, 2.2, 3.2};
    for (int i=0; i<3; ++i) {
      B.getVectorNonConst(i)->putScalar(bc[i]);
    }

    Eigen::MatrixXd C(A.getNumVectors(), B.getNumVectors());
    C.setConstant(NAN);

    // C = 1*C + 1.5 A^T B
    pressio::ops::product(
        pressio::transpose(),
        pressio::nontranspose(),
        1.5, A, B, 0.0, C);

    for (auto i=0; i<C.rows(); i++){
        for (auto j=0; j<C.cols(); j++){
            const auto gold = ac[i]*A.getGlobalLength()*1.5*bc[j];
            EXPECT_NEAR( C(i,j), gold, 1e-12);
        }
    }
}

TEST_F(tpetraMultiVectorGlobSize15Fixture, mv_T_self_storein_eigen_C)
{
    auto A = pressio::ops::clone(*myMv_);
    std::array<double, 4> ac{1.,2.,3.,4.};
    for (std::size_t i=0; i<A.getNumVectors(); ++i) {
      A.getVectorNonConst(i)->putScalar(ac[i]);
    }

    Eigen::MatrixXd C(A.getNumVectors(), A.getNumVectors());
    C.setConstant(0.);

    // C = 1*C + 1.5 A^T A
    pressio::ops::product(
        pressio::transpose(),
        pressio::nontranspose(),
        1.5, A, 1.0, C);

    if(rank_==0){
        std::cout << C << std::endl;
    }

    for (auto i=0; i<C.rows(); i++){
        for (auto j=0; j<C.cols(); j++){
            const auto gold = ac[i]*A.getGlobalLength()*1.5*ac[j];
            EXPECT_NEAR( C(i,j), gold, 1e-12);
        }
    }
}

TEST_F(tpetraMultiVectorGlobSize15Fixture, mv_T_self_storein_eigen_C_beta0)
{
    auto A = pressio::ops::clone(*myMv_);
    std::array<double, 4> ac{1.,2.,3.,4.};
    for (std::size_t i=0; i<A.getNumVectors(); ++i) {
      A.getVectorNonConst(i)->putScalar(ac[i]);
    }

    Eigen::MatrixXd C(A.getNumVectors(), A.getNumVectors());
    C.setConstant(NAN);

    // C = 0*NAN + 1.5 A^T A
    pressio::ops::product(
        pressio::transpose(),
        pressio::nontranspose(),
        1.5, A, 0.0, C);

    if(rank_==0){
        std::cout << C << std::endl;
    }

    for (auto i=0; i<C.rows(); i++){
        for (auto j=0; j<C.cols(); j++){
            const auto gold = ac[i]*A.getGlobalLength()*1.5*ac[j];
            EXPECT_NEAR( C(i,j), gold, 1e-12);
        }
    }
}

TEST_F(tpetraMultiVectorGlobSize15Fixture, mv_T_self_storein_eigen_C_test2)
{
  /*
    C = 2 * C + 3.*A^T A

    A[:,0] = 1.
    A[:,1] = 2.
    A[:,2] = 3.
    A[:,3] = 4.

    C = [1  2   3  4;
         5  6   7  8;
	 9 10  11  12;
	 13 14 15  16]
   */

    auto A = pressio::ops::clone(*myMv_);
    std::array<double, 4> ac{1.,2.,3.,4.};
    for (std::size_t i=0; i<A.getNumVectors(); ++i) {
      A.getVectorNonConst(i)->putScalar(ac[i]);
    }

    Eigen::MatrixXd C(A.getNumVectors(), A.getNumVectors());
    double s=1.;
    for (auto i=0; i<C.rows(); i++){
        for (auto j=0; j<C.cols(); j++){
	  C(i,j) = s;
	  s += 1.;
	}
    }

    pressio::ops::product(
        pressio::transpose(),
        pressio::nontranspose(),
        2., A, 3.0, C);

    if(rank_==0){
        std::cout << C << std::endl;
    }

    EXPECT_NEAR( C(0,0), 3.*1. + 2.*15., 1e-12);
    EXPECT_NEAR( C(0,1), 3.*2. + 2.*30., 1e-12);
    EXPECT_NEAR( C(0,2), 3.*3. + 2.*45., 1e-12);
    EXPECT_NEAR( C(0,3), 3.*4. + 2.*60., 1e-12);

    EXPECT_NEAR( C(1,0), 3.*5. + 2.*30., 1e-12);
    EXPECT_NEAR( C(1,1), 3.*6. + 2.*60., 1e-12);
    EXPECT_NEAR( C(1,2), 3.*7. + 2.*90., 1e-12);
    EXPECT_NEAR( C(1,3), 3.*8. + 2.*120., 1e-12);

    EXPECT_NEAR( C(2,0), 3.*9.  + 2.*45., 1e-12);
    EXPECT_NEAR( C(2,1), 3.*10. + 2.*90., 1e-12);
    EXPECT_NEAR( C(2,2), 3.*11. + 2.*135., 1e-12);
    EXPECT_NEAR( C(2,3), 3.*12. + 2.*180., 1e-12);

    EXPECT_NEAR( C(3,0), 3.*13. + 2.*60., 1e-12);
    EXPECT_NEAR( C(3,1), 3.*14. + 2.*120., 1e-12);
    EXPECT_NEAR( C(3,2), 3.*15. + 2.*180., 1e-12);
    EXPECT_NEAR( C(3,3), 3.*16. + 2.*240., 1e-12);
}
#endif


TEST_F(tpetraMultiVectorGlobSize15Fixture, mv_T_self_storein_kokkos_C)
{
    auto A = pressio::ops::clone(*myMv_);
    std::array<double, 4> ac{1.,2.,3.,4.};
    for (std::size_t i=0; i<A.getNumVectors(); ++i) {
      A.getVectorNonConst(i)->putScalar(ac[i]);
    }

    Kokkos::View<double**, Kokkos::LayoutLeft> C("C", A.getNumVectors(), A.getNumVectors());

    // C = 1*C + 1.5 A^T A
    pressio::ops::product(
        pressio::transpose(),
        pressio::nontranspose(),
        1.5, A, 1.0, C);

    auto C2 = pressio::ops::product<Kokkos::View<double**, Kokkos::LayoutLeft>>(
        pressio::transpose(), pressio::nontranspose(), 1.5, A);

    auto C_h = Kokkos::create_mirror_view(C);
    auto C2_h = Kokkos::create_mirror_view(C2);
    for (std::size_t i=0; i<C.extent(0); i++){
        for (std::size_t j=0; j<C.extent(1); j++){
            const auto gold = ac[i]*A.getGlobalLength()*1.5*ac[j];
            EXPECT_NEAR( C_h(i,j), gold, 1e-12);
            EXPECT_NEAR( C2_h(i,j), gold, 1e-12);
        }
    }
}

TEST_F(tpetraMultiVectorGlobSize15Fixture, mv_T_self_storein_kokkos_C_beta0)
{
    auto A = pressio::ops::clone(*myMv_);
    std::array<double, 4> ac{1.,2.,3.,4.};
    for (std::size_t i=0; i<A.getNumVectors(); ++i) {
      A.getVectorNonConst(i)->putScalar(ac[i]);
    }

    Kokkos::View<double**, Kokkos::LayoutLeft> C("C", A.getNumVectors(), A.getNumVectors());
    Kokkos::deep_copy(C, NAN);

    // C = 0*NAN + 1.5 A^T A
    pressio::ops::product(
        pressio::transpose(),
        pressio::nontranspose(),
        1.5, A, 0.0, C);

    auto C2 = pressio::ops::product<Kokkos::View<double**, Kokkos::LayoutLeft>>(
        pressio::transpose(), pressio::nontranspose(), 1.5, A);

    auto C_h = Kokkos::create_mirror_view(C);
    auto C2_h = Kokkos::create_mirror_view(C2);
    for (std::size_t i=0; i<C.extent(0); i++){
        for (std::size_t j=0; j<C.extent(1); j++){
            const auto gold = ac[i]*A.getGlobalLength()*1.5*ac[j];
            EXPECT_NEAR( C_h(i,j), gold, 1e-12);
            EXPECT_NEAR( C2_h(i,j), gold, 1e-12);
        }
    }
}


TEST_F(tpetraMultiVectorGlobSize15Fixture, mv_T_mv_storein_kokkos_C)
{
    auto A = pressio::ops::clone(*myMv_);
    std::array<double, 4> ac{1.,2.,3.,4.};
    for (std::size_t i=0; i<A.getNumVectors(); ++i) {
      A.getVectorNonConst(i)->putScalar(ac[i]);
    }

    mvec_t B(contigMap_, 3);
    std::array<double, 3> bc{1.2, 2.2, 3.2};
    for (int i=0; i<3; ++i) {
      B.getVectorNonConst(i)->putScalar(bc[i]);
    }

    Kokkos::View<double**, Kokkos::LayoutLeft> C("C", A.getNumVectors(), B.getNumVectors());

    // C = 1*C + 1.5 A^T B
    pressio::ops::product(
        pressio::transpose(),
        pressio::nontranspose(),
        1.5, A, B, 1.0, C);

    auto C_h = Kokkos::create_mirror_view(C);
    for (std::size_t i=0; i<C.extent(0); i++){
        for (std::size_t j=0; j<C.extent(1); j++){
            const auto gold = ac[i]*A.getGlobalLength()*1.5*bc[j];
            EXPECT_NEAR( C_h(i,j), gold, 1e-12);
        }
    }
}

TEST_F(tpetraMultiVectorGlobSize15Fixture, mv_T_mv_storein_kokkos_C_beta0)
{
    auto A = pressio::ops::clone(*myMv_);
    std::array<double, 4> ac{1.,2.,3.,4.};
    for (std::size_t i=0; i<A.getNumVectors(); ++i) {
      A.getVectorNonConst(i)->putScalar(ac[i]);
    }

    mvec_t B(contigMap_, 3);
    std::array<double, 3> bc{1.2, 2.2, 3.2};
    for (int i=0; i<3; ++i) {
      B.getVectorNonConst(i)->putScalar(bc[i]);
    }

    Kokkos::View<double**, Kokkos::LayoutLeft> C("C", A.getNumVectors(), B.getNumVectors());
    Kokkos::deep_copy(C, NAN);

    // C = 0*NAN + 1.5 A^T B
    pressio::ops::product(
        pressio::transpose(),
        pressio::nontranspose(),
        1.5, A, B, 0.0, C);

    auto C_h = Kokkos::create_mirror_view(C);
    for (std::size_t i=0; i<C.extent(0); i++){
        for (std::size_t j=0; j<C.extent(1); j++){
            const auto gold = ac[i]*A.getGlobalLength()*1.5*bc[j];
            EXPECT_NEAR( C_h(i,j), gold, 1e-12);
        }
    }
}
