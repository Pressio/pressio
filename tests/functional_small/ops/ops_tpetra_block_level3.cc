
#include "tpetra_block_only_fixtures.hpp"
#include "pressio/ops.hpp"

#ifdef PRESSIO_ENABLE_TPL_EIGEN
TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       mv_T_mv_storein_eigen_C)
{
  auto A = pressio::ops::clone(*myMv_);
  std::array<double, 4> ac{1.,2.,3.,4.};
  for (std::size_t i=0; i<A.getMultiVectorView().getNumVectors(); ++i) {
    A.getMultiVectorView().getVectorNonConst(i)->putScalar(ac[i]);
  }

  mvec_t B(*contigMap_, blockSize_, 4);
  std::array<double, 4> bc{1.2, 2.2, 3.2, -4.1};
  for (int i=0; i<4; ++i) {
    B.getMultiVectorView().getVectorNonConst(i)->putScalar(bc[i]);
  }

  Eigen::MatrixXd C(A.getNumVectors(), B.getNumVectors());
  C.setConstant(0.);

  // C = 1*C + 1.5 A^T B
  pressio::ops::product(pressio::transpose(),
			pressio::nontranspose(),
			1.5, A, B, 1.0, C);

  for (auto i=0; i<C.rows(); i++){
    for (auto j=0; j<C.cols(); j++){
      const auto gold = ac[i]*A.getMultiVectorView().getGlobalLength()*1.5*bc[j];
      EXPECT_NEAR( C(i,j), gold, 1e-12);
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       mv_T_self_storein_eigen_C)
{
  auto A = pressio::ops::clone(*myMv_);
  std::array<double, 4> ac{1.,2.,3.,4.};
  for (std::size_t i=0; i<A.getMultiVectorView().getNumVectors(); ++i) {
    A.getMultiVectorView().getVectorNonConst(i)->putScalar(ac[i]);
  }

  Eigen::MatrixXd C(A.getNumVectors(), A.getNumVectors());
  C.setConstant(0.);

  // C = 1*C + 1.5 A^T A
  pressio::ops::product(pressio::transpose(),
			pressio::nontranspose(),
			1.5, A, 1.0, C);

  if(rank_==0){
    std::cout << C << std::endl;
  }

  for (auto i=0; i<C.rows(); i++){
    for (auto j=0; j<C.cols(); j++){
      const auto gold = ac[i]*A.getMultiVectorView().getGlobalLength()*1.5*ac[j];
      EXPECT_NEAR( C(i,j), gold, 1e-12);
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       mv_T_self_create_result_eigen_C)
{
  auto A = pressio::ops::clone(*myMv_);
  std::array<double, 4> ac{1.,2.,3.,4.};
  for (std::size_t i=0; i<A.getMultiVectorView().getNumVectors(); ++i) {
    A.getMultiVectorView().getVectorNonConst(i)->putScalar(ac[i]);
  }

  // C = 1.5 A^T A
  auto C = pressio::ops::product<Eigen::MatrixXd>(pressio::transpose(),
						  pressio::nontranspose(),
						  1.5, A);

  if(rank_==0){
    std::cout << C << std::endl;
  }

  for (auto i=0; i<C.rows(); i++){
    for (auto j=0; j<C.cols(); j++){
      const auto gold = ac[i]*A.getMultiVectorView().getGlobalLength()*1.5*ac[j];
      EXPECT_NEAR( C(i,j), gold, 1e-12);
    }
  }
}
#endif


TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       mv_T_self_storein_kokkos_C)
{
  auto A = pressio::ops::clone(*myMv_);
  std::array<double, 4> ac{1.,2.,3.,4.};
  for (std::size_t i=0; i<A.getNumVectors(); ++i) {
    A.getMultiVectorView().getVectorNonConst(i)->putScalar(ac[i]);
  }

  Kokkos::View<double**, Kokkos::LayoutLeft> C("C", A.getNumVectors(), A.getNumVectors());

  // C = 1*C + 1.5 A^T A
  pressio::ops::product(pressio::transpose(),
			pressio::nontranspose(),
			1.5, A, 1.0, C);

  auto C2 = pressio::ops::product<
    Kokkos::View<double**, Kokkos::LayoutLeft>>(pressio::transpose(),
						pressio::nontranspose(), 1.5, A);

  auto C_h = Kokkos::create_mirror_view(C);
  auto C2_h = Kokkos::create_mirror_view(C2);
  for (std::size_t i=0; i<C.extent(0); i++){
    for (std::size_t j=0; j<C.extent(1); j++){
      const auto gold = ac[i]*A.getMultiVectorView().getGlobalLength()*1.5*ac[j];
      EXPECT_NEAR( C_h(i,j), gold, 1e-12);
      EXPECT_NEAR( C2_h(i,j), gold, 1e-12);
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       mv_T_mv_storein_kokkos_C)
{
  auto A = pressio::ops::clone(*myMv_);
  std::array<double, 4> ac{1.,2.,3.,4.};
  for (std::size_t i=0; i<A.getNumVectors(); ++i) {
    A.getMultiVectorView().getVectorNonConst(i)->putScalar(ac[i]);
  }

  mvec_t B(*contigMap_, blockSize_, 3);
  std::array<double, 3> bc{1.2, 2.2, 3.2};
  for (int i=0; i<3; ++i) {
    B.getMultiVectorView().getVectorNonConst(i)->putScalar(bc[i]);
  }

  Kokkos::View<double**, Kokkos::LayoutLeft> C("C", A.getNumVectors(), B.getNumVectors());

  // C = 1*C + 1.5 A^T B
  pressio::ops::product(pressio::transpose(),
			pressio::nontranspose(),
			1.5, A, B, 1.0, C);

  auto C_h = Kokkos::create_mirror_view(C);
  for (std::size_t i=0; i<C.extent(0); i++){
    for (std::size_t j=0; j<C.extent(1); j++){
      const auto gold = ac[i]*A.getMultiVectorView().getGlobalLength()*1.5*bc[j];
      EXPECT_NEAR( C_h(i,j), gold, 1e-12);
    }
  }
}
