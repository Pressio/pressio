
#include "tpetra_block_only_fixtures.hpp"
#include "pressio/ops.hpp"

#ifdef PRESSIO_ENABLE_TPL_EIGEN
TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       mv_T_mv_storein_eigen_C)
{
  // C = 1*C + 1.5 A^T B
  pressio::ops::product(pressio::transpose(),
			pressio::nontranspose(),
			1.5, A, *B, 1.0, C_eigen);

  for (auto i = 0; i < C_eigen.rows(); i++){
    for (auto j = 0; j < C_eigen.cols(); j++){
      EXPECT_NEAR(C_eigen(i, j), gold_ab(i, j), 1e-12);
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       mv_T_self_storein_eigen_C)
{
  // C = 1*C + 1.5 A^T A
  pressio::ops::product(pressio::transpose(),
			pressio::nontranspose(),
			1.5, A, 1.0, C_eigen);

  if(rank_==0){
    std::cout << C_eigen << std::endl;
  }

  for (auto i = 0; i < C_eigen.rows(); i++){
    for (auto j = 0; j < C_eigen.cols(); j++){
      EXPECT_NEAR(C_eigen(i, j), gold_a(i, j), 1e-12);
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       mv_T_self_create_result_eigen_C)
{
  // C = 1.5 A^T A
  auto C = pressio::ops::product<Eigen::MatrixXd>(pressio::transpose(),
						  pressio::nontranspose(),
						  1.5, A);

  auto C2 = pressio::ops::product<Eigen::MatrixXd>(pressio::transpose(),
						  pressio::nontranspose(),
						  1.5, A, *B);

  if(rank_==0){
    std::cout << C << "\n\n"
              << C2 << std::endl;
  }

  for (auto i=0; i<C.rows(); i++){
    for (auto j=0; j<C.cols(); j++){
      EXPECT_NEAR(C(i, j), gold_a(i, j), 1e-12);
      EXPECT_NEAR(C2(i, j), gold_ab(i, j), 1e-12);
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       mv_T_self_storein_eigen_C_beta0)
{
  C_eigen.setConstant(std::nan("0"));

  // C = 0*NAN + 1.5 A^T A
  pressio::ops::product(pressio::transpose(),
			pressio::nontranspose(),
			1.5, A, 0.0, C_eigen);

  if(rank_==0){
    std::cout << C_eigen << std::endl;
  }

  for (auto i = 0; i < C_eigen.rows(); i++){
    for (auto j = 0; j < C_eigen.cols(); j++){
      EXPECT_NEAR(C_eigen(i, j), gold_a(i, j), 1e-12);
    }
  }
}
#endif

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       mv_T_self_storein_kokkos_C)
{
  // C = 1*C + 1.5 A^T A
  pressio::ops::product(pressio::transpose(),
			pressio::nontranspose(),
			1.5, A, 1.0, C_kokkos);

  auto C2_kokkos = pressio::ops::product<
    Kokkos::View<double**, Kokkos::LayoutLeft>>(pressio::transpose(),
						pressio::nontranspose(), 1.5, A);

  auto C_h = Kokkos::create_mirror_view(C_kokkos);
  auto C2_h = Kokkos::create_mirror_view(C2_kokkos);
  for (std::size_t i = 0; i < C_kokkos.extent(0); i++){
    for (std::size_t j = 0; j < C_kokkos.extent(1); j++){
      EXPECT_NEAR(C_h(i, j), gold_a(i, j), 1e-12);
      EXPECT_NEAR(C2_h(i, j), gold_a(i, j), 1e-12);
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       mv_T_self_storein_kokkos_C_beta0)
{
  Kokkos::deep_copy(C_kokkos, std::nan("0"));

  // C = 0*NAN + 1.5 A^T A
  pressio::ops::product(pressio::transpose(),
			pressio::nontranspose(),
			1.5, A, 0.0, C_kokkos);

  auto C2_kokkos = pressio::ops::product<
    Kokkos::View<double**, Kokkos::LayoutLeft>>(pressio::transpose(),
						pressio::nontranspose(), 1.5, A);

  auto C_h = Kokkos::create_mirror_view(C_kokkos);
  auto C2_h = Kokkos::create_mirror_view(C2_kokkos);
  for (std::size_t i = 0; i < C_kokkos.extent(0); i++){
    for (std::size_t j = 0; j < C_kokkos.extent(1); j++){
      EXPECT_NEAR(C_h(i, j), gold_a(i, j), 1e-12);
      EXPECT_NEAR(C2_h(i, j), gold_a(i, j), 1e-12);
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       mv_T_mv_storein_kokkos_C)
{
  // C = 1*C + 1.5 A^T B
  pressio::ops::product(pressio::transpose(),
			pressio::nontranspose(),
			1.5, A, *B, 1.0, C_kokkos);

  auto C_h = Kokkos::create_mirror_view(C_kokkos);
  for (std::size_t i = 0; i < C_kokkos.extent(0); i++){
    for (std::size_t j = 0; j < C_kokkos.extent(1); j++){
      EXPECT_NEAR(C_h(i, j), gold_ab(i, j), 1e-12);
    }
  }
}

TEST_F(tpetraBlockMultiVectorGlobSize15NVec3BlockSize4Fixture,
       mv_T_mv_storein_kokkos_C_beta0)
{
  Kokkos::deep_copy(C_kokkos, std::nan("0"));

  // C = 0*NAN + 1.5 A^T B
  pressio::ops::product(pressio::transpose(),
			pressio::nontranspose(),
			1.5, A, *B, 0.0, C_kokkos);

  auto C_h = Kokkos::create_mirror_view(C_kokkos);
  for (std::size_t i = 0; i < C_kokkos.extent(0); i++){
    for (std::size_t j = 0; j < C_kokkos.extent(1); j++){
      EXPECT_NEAR(C_h(i, j), gold_ab(i, j), 1e-12);
    }
  }
}
