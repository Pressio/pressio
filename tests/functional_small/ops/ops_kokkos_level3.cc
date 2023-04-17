
#include <gtest/gtest.h>
#include "pressio/ops.hpp"

using mat_t = Kokkos::View<double**>;

#define OPS_KOKKOS_DENSE_MAT_MAT_PROD(OPERAND) \
  mat_t M("M",4,3); \
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
  M_h(3,1) = 2.; \
  M_h(3,2) = 2.; \
  Kokkos::deep_copy(M, M_h); \
  mat_t myR("MYR", 4,4);         \
  myR(0, 0) = std::nan("0"); /* simulate uninitialized NaN */ \
  constexpr auto beta  = ::pressio::utils::Constants<double>::zero(); \
  constexpr auto alpha = ::pressio::utils::Constants<double>::one();  \
  pressio::ops::product(pressio::nontranspose(), pressio::nontranspose(), \
      alpha, M, OPERAND, beta, myR);\
  auto myR_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), myR); \
  EXPECT_DOUBLE_EQ(myR_h(0,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR_h(0,1), 3.0); \
  EXPECT_DOUBLE_EQ(myR_h(0,2), 6.0); \
  EXPECT_DOUBLE_EQ(myR_h(0,3), 9.0); \
  EXPECT_DOUBLE_EQ(myR_h(1,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR_h(1,1), 6.0); \
  EXPECT_DOUBLE_EQ(myR_h(1,2), 12.0); \
  EXPECT_DOUBLE_EQ(myR_h(1,3), 18.0); \
  EXPECT_DOUBLE_EQ(myR_h(2,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR_h(2,1), 1.0); \
  EXPECT_DOUBLE_EQ(myR_h(2,2), 2.0); \
  EXPECT_DOUBLE_EQ(myR_h(2,3), 3.0); \
  EXPECT_DOUBLE_EQ(myR_h(3,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR_h(3,1), 6.0); \
  EXPECT_DOUBLE_EQ(myR_h(3,2), 12.0); \
  EXPECT_DOUBLE_EQ(myR_h(3,3), 18.0); \


#define OPS_KOKKOS_DENSE_MAT_T_MAT_PROD(OPERAND) \
  mat_t M("M",4,3); \
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
  M_h(3,1) = 2.; \
  M_h(3,2) = 2.; \
  Kokkos::deep_copy(M, M_h); \
  mat_t myR("MYR", 3,3);         \
  myR(0, 0) = std::nan("0"); /* simulate uninitialized NaN */ \
  constexpr auto beta  = ::pressio::utils::Constants<double>::zero(); \
  constexpr auto alpha = ::pressio::utils::Constants<double>::one();  \
  pressio::ops::product(pressio::transpose(), pressio::nontranspose(), \
      alpha, M, OPERAND, beta, myR); \
  auto myR_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), myR); \
  EXPECT_DOUBLE_EQ(myR_h(0,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR_h(0,1), 5.0); \
  EXPECT_DOUBLE_EQ(myR_h(0,2), 10.0); \
  EXPECT_DOUBLE_EQ(myR_h(1,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR_h(1,1), 3.0); \
  EXPECT_DOUBLE_EQ(myR_h(1,2), 6.0); \
  EXPECT_DOUBLE_EQ(myR_h(2,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR_h(2,1), 8.0); \
  EXPECT_DOUBLE_EQ(myR_h(2,2), 16.0); \


#define OPS_KOKKOS_DENSE_MAT_T_SELF_PROD \
  mat_t M("M",4,3); \
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
  M_h(3,1) = 2.; \
  M_h(3,2) = 2.; \
  Kokkos::deep_copy(M, M_h); \
  mat_t myR("MYR", 3,3);         \
  myR(0, 0) = std::nan("0"); /* simulate uninitialized NaN */ \
  constexpr auto beta  = ::pressio::utils::Constants<double>::zero(); \
  constexpr auto alpha = ::pressio::utils::Constants<double>::one();  \
  pressio::ops::product(pressio::transpose(), pressio::nontranspose(), alpha, M, beta, myR); \
  auto myR_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), myR); \
  EXPECT_DOUBLE_EQ(myR_h(0,0), 9.0); \
  EXPECT_DOUBLE_EQ(myR_h(0,1), 6.0); \
  EXPECT_DOUBLE_EQ(myR_h(0,2), 12.0); \
  EXPECT_DOUBLE_EQ(myR_h(1,0), 6.0); \
  EXPECT_DOUBLE_EQ(myR_h(1,1), 5.0); \
  EXPECT_DOUBLE_EQ(myR_h(1,2), 7.0); \
  EXPECT_DOUBLE_EQ(myR_h(2,0), 12.0); \
  EXPECT_DOUBLE_EQ(myR_h(2,1), 7.0); \
  EXPECT_DOUBLE_EQ(myR_h(2,2), 18.0); \
  auto myR1 = pressio::ops::product<mat_t>(pressio::transpose(), pressio::nontranspose(), alpha, M); \
  auto myR1_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), myR1); \
  EXPECT_DOUBLE_EQ(myR1_h(0,0), 9.0); \
  EXPECT_DOUBLE_EQ(myR1_h(0,1), 6.0); \
  EXPECT_DOUBLE_EQ(myR1_h(0,2), 12.0); \
  EXPECT_DOUBLE_EQ(myR1_h(1,0), 6.0); \
  EXPECT_DOUBLE_EQ(myR1_h(1,1), 5.0); \
  EXPECT_DOUBLE_EQ(myR1_h(1,2), 7.0); \
  EXPECT_DOUBLE_EQ(myR1_h(2,0), 12.0); \
  EXPECT_DOUBLE_EQ(myR1_h(2,1), 7.0); \
  EXPECT_DOUBLE_EQ(myR1_h(2,2), 18.0); \

namespace {
template<class T>
void fillOperand1(T & M)
{
  auto M_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), M);
  for (int i=0; i<3; ++i){
    for (int j=0; j<4; ++j){
      M_h(i,j) = (double) j;
    }
  }
  Kokkos::deep_copy(M, M_h);
}
template<class T>
void fillOperand2(T & M)
{
  auto M_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), M);
  for (int i=0; i<4; ++i){
    for (int j=0; j<3; ++j){
      M_h(i,j) = (double) j;
    }
  }
  Kokkos::deep_copy(M, M_h);
}
}//end namespace

TEST(ops_kokkos, dense_matrix_dense_matrix_prod)
{
  mat_t operand("M", 3,4);
  fillOperand1(operand);
  OPS_KOKKOS_DENSE_MAT_MAT_PROD(operand);
}

TEST(ops_kokkos, dense_matrix_T_dense_matrix_prod)
{
  mat_t operand("M", 4,3);
  fillOperand2(operand);
  OPS_KOKKOS_DENSE_MAT_T_MAT_PROD(operand);
}

TEST(ops_kokkos, dense_matrix_T_self_prod)
{
  OPS_KOKKOS_DENSE_MAT_T_SELF_PROD;
}

