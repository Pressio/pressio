
#include <gtest/gtest.h>
#include "pressio/ops.hpp"

#define OPS_KOKKOS_DENSEMATRIX_VEC_PROD(VECIN) \
  Kokkos::View<double**> M("M",3,3); \
  auto M_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), M); \
  M_h(0,0) = 1.; \
  M_h(0,1) = 0.; \
  M_h(0,2) = 2.; \
  M_h(1,0) = 2.; \
  M_h(1,1) = 1.; \
  M_h(1,2) = 3.; \
  M_h(2,0) = 0.; \
  M_h(2,1) = 0.; \
  M_h(2,2) = 1.; \
  Kokkos::deep_copy(M, M_h); \
  Kokkos::View<double*> myR("myR", 3); \
  myR(0) = NAN; /* simulate uninitialized Nan */ \
  constexpr auto beta  = ::pressio::utils::Constants<double>::zero(); \
  constexpr auto alpha = ::pressio::utils::Constants<double>::one();  \
  pressio::ops::product(pressio::nontranspose(), alpha, M, VECIN, beta, myR); \
  auto myR_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), myR); \
  EXPECT_DOUBLE_EQ( myR_h(0), 16.0); \
  EXPECT_DOUBLE_EQ( myR_h(1), 28.0); \
  EXPECT_DOUBLE_EQ( myR_h(2), 6.0);  \

#define OPS_KOKKOS_DENSEMATRIX_T_VEC_PROD(VECIN) \
  Kokkos::View<double**> M("M",4,3); \
  auto M_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), M); \
  M_h(0,0) = 1.; \
  M_h(0,1) = 0.; \
  M_h(0,2) = 2.; \
  M_h(1,0) = 2.; \
  M_h(1,1) = 1.; \
  M_h(1,2) = 3.; \
  M_h(2,0) = 0.; \
  M_h(2,1) = 0.; \
  M_h(2,2) = 1.; \
  M_h(3,0) = 2.; \
  M_h(3,1) = 3.; \
  M_h(3,2) = 4.; \
  Kokkos::deep_copy(M, M_h); \
  Kokkos::View<double*> myR("myR", 3); \
  myR(0) = NAN; /* simulate uninitialized Nan */ \
  constexpr auto beta  = ::pressio::utils::Constants<double>::zero(); \
  constexpr auto alpha = ::pressio::utils::Constants<double>::one();  \
  pressio::ops::product(pressio::transpose(), alpha, M, VECIN, beta, myR); \
  auto myR_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), myR); \
  EXPECT_DOUBLE_EQ( myR_h(0), 14.);  \
  EXPECT_DOUBLE_EQ( myR_h(1), 11.0); \
  EXPECT_DOUBLE_EQ( myR_h(2), 8.+6.+6.+12.);  \


TEST(ops_kokkos, dense_matrix_vector_prod)
{
  Kokkos::View<double*> a("a", 3);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  a_h(0) = 4.; a_h(1) = 2.; a_h(2) = 6.;
  Kokkos::deep_copy(a, a_h);
  OPS_KOKKOS_DENSEMATRIX_VEC_PROD(a);
}

TEST(ops_kokkos, dense_matrix_T_vector_prod)
{
  Kokkos::View<double*> a("a", 4);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  a_h(0) = 4.; a_h(1) = 2.; a_h(2) = 6.; a_h(3)=3.;
  Kokkos::deep_copy(a, a_h);
  OPS_KOKKOS_DENSEMATRIX_T_VEC_PROD(a);
}

TEST(ops_kokkos, dense_matrix_span_prod)
{
  Kokkos::View<double*> a("a", 7);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  a_h(3)=4.;
  a_h(4)=2.;
  a_h(5)=6.;
  Kokkos::deep_copy(a, a_h);

  const auto sp = pressio::span(a, 3, 3);
  OPS_KOKKOS_DENSEMATRIX_VEC_PROD(sp);
}

TEST(ops_kokkos, dense_matrix_T_span_prod)
{
  Kokkos::View<double*> a("a", 7);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  a_h(3)=4.;
  a_h(4)=2.;
  a_h(5)=6.;
  a_h(6)=3.;
  Kokkos::deep_copy(a, a_h);

  const auto sp = pressio::span(a, 3, 4);
  OPS_KOKKOS_DENSEMATRIX_T_VEC_PROD(sp);
}

TEST(ops_kokkos, dense_matrix_diag_prod)
{
  Kokkos::View<double**> M0("M0", 3,3);
  auto M0_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), M0);
  M0_h(0,0)=4.;
  M0_h(1,1)=2.;
  M0_h(2,2)=6.;
  Kokkos::deep_copy(M0, M0_h);

  const auto exp = pressio::diag(M0);
  OPS_KOKKOS_DENSEMATRIX_VEC_PROD(exp);
}

TEST(ops_kokkos, dense_matrix_T_diag_prod)
{
  Kokkos::View<double**> M0("M0", 4,4);
  auto M0_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), M0);
  M0_h(0,0)=4.;
  M0_h(1,1)=2.;
  M0_h(2,2)=6.;
  M0_h(3,3)=3.;
  Kokkos::deep_copy(M0, M0_h);

  const auto exp = pressio::diag(M0);
  OPS_KOKKOS_DENSEMATRIX_T_VEC_PROD(exp);
}

