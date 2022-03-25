
#include <gtest/gtest.h>
#include "pressio/ops.hpp"

#define OPS_EIGEN_DENSEMATRIX_VEC_PROD(VECIN) \
  using M_t = Eigen::MatrixXd; \
  M_t M(3,3);                  \
  M << 1.,0.,2.,2.,1.,3.,0.,0.,1.;      \
  Eigen::VectorXd myR(3);      \
  myR(0) = NAN; /* simulate uninitialized Nan */ \
  constexpr auto beta  = ::pressio::utils::Constants<double>::zero(); \
  constexpr auto alpha = ::pressio::utils::Constants<double>::one();  \
  pressio::ops::product(pressio::nontranspose(), alpha, M, VECIN, beta, myR); \
  EXPECT_DOUBLE_EQ( myR(0), 16.0); \
  EXPECT_DOUBLE_EQ( myR(1), 28.0); \
  EXPECT_DOUBLE_EQ( myR(2), 6.0);  \
  /* also cover non-zero beta */ \
  constexpr auto beta1 = ::pressio::utils::Constants<double>::one(); \
  pressio::ops::product(pressio::nontranspose(), alpha, M, VECIN, beta1, myR); \
  EXPECT_DOUBLE_EQ( myR(0), 32.0); \
  EXPECT_DOUBLE_EQ( myR(1), 56.0); \
  EXPECT_DOUBLE_EQ( myR(2), 12.0);  \

#define OPS_EIGEN_DENSEMATRIX_T_VEC_PROD(VECIN) \
  using M_t = Eigen::MatrixXd; \
  M_t M(4,3);                  \
  M << 1.,0.,2.,2.,1.,3.,0.,0.,1.,2.,3.,4.; \
                               \
  Eigen::VectorXd myR(3);      \
  myR(0) = NAN; /* simulate uninitialized Nan */ \
  constexpr auto beta  = ::pressio::utils::Constants<double>::zero(); \
  constexpr auto alpha = ::pressio::utils::Constants<double>::one();  \
  pressio::ops::product(pressio::transpose(), alpha, M, VECIN, beta, myR); \
  EXPECT_DOUBLE_EQ( myR(0), 14.);  \
  EXPECT_DOUBLE_EQ( myR(1), 11.0); \
  EXPECT_DOUBLE_EQ( myR(2), 8.+6.+6.+12.);  \
  /* also cover non-zero beta */ \
  constexpr auto beta1 = ::pressio::utils::Constants<double>::one(); \
  pressio::ops::product(pressio::transpose(), alpha, M, VECIN, beta1, myR); \
  EXPECT_DOUBLE_EQ( myR(0), 28.0); \
  EXPECT_DOUBLE_EQ( myR(1), 22.0); \
  EXPECT_DOUBLE_EQ( myR(2), 64.0);  \
  auto y = pressio::ops::product<Eigen::VectorXd>(pressio::transpose(), alpha, M, VECIN); \
  EXPECT_DOUBLE_EQ( y(0), 14.);  \
  EXPECT_DOUBLE_EQ( y(1), 11.0); \
  EXPECT_DOUBLE_EQ( y(2), 8.+6.+6.+12.);  \


TEST(ops_eigen, dense_matrix_vector_prod)
{
  using V_t = Eigen::VectorXd;
  V_t a(3); a << 4.,2.,6;
  OPS_EIGEN_DENSEMATRIX_VEC_PROD(a);
}

TEST(ops_eigen, dense_matrix_T_vector_prod)
{
  using V_t = Eigen::VectorXd;
  V_t a(4); a << 4.,2.,6,3.;
  OPS_EIGEN_DENSEMATRIX_T_VEC_PROD(a);
}

TEST(ops_eigen, dense_matrix_span_prod)
{
  using V_t = Eigen::VectorXd;
  V_t a(7); 
  a(3)=4.;
  a(4)=2.;
  a(5)=6.;

  const auto sp = pressio::span(a, 3, 3);
  OPS_EIGEN_DENSEMATRIX_VEC_PROD(sp);
}


TEST(ops_eigen, dense_matrix_T_span_prod)
{
  using V_t = Eigen::VectorXd;
  V_t a(7);
  a(3)=4.;
  a(4)=2.;
  a(5)=6.;
  a(6)=3.;

  const auto sp = pressio::span(a, 3, 4);
  OPS_EIGEN_DENSEMATRIX_T_VEC_PROD(sp);
}

TEST(ops_eigen, dense_matrix_diag_prod)
{
  using T = Eigen::MatrixXd;
  T M0(3,3); 
  M0(0,0)=4.;
  M0(1,1)=2.;
  M0(2,2)=6.;

  const auto exp = pressio::diag(M0);
  OPS_EIGEN_DENSEMATRIX_VEC_PROD(exp);
}

TEST(ops_eigen, dense_matrix_T_diag_prod)
{
  using T = Eigen::MatrixXd;
  T M0(4,4);
  M0(0,0)=4.;
  M0(1,1)=2.;
  M0(2,2)=6.;
  M0(3,3)=3.;

  const auto exp = pressio::diag(M0);
  OPS_EIGEN_DENSEMATRIX_T_VEC_PROD(exp);
}

