
#include <gtest/gtest.h>
#include "pressio/ops.hpp"

#define OPS_EIGEN_DENSEMATRIX_VEC_PROD(VECIN) \
  using M_t = Eigen::MatrixXd; \
  M_t M(3,3);                  \
  M << 1.,0.,2.,2.,1.,3.,0.,0.,1.;      \
  Eigen::VectorXd myR(3);      \
  myR(0) = std::nan("0"); /* simulate uninitialized NaN */ \
  constexpr auto beta0  = static_cast<double>(0); \
  constexpr auto alpha1 = static_cast<double>(1);  \
  pressio::ops::product(pressio::nontranspose(), alpha1, M, VECIN, beta0, myR); \
  EXPECT_DOUBLE_EQ( myR(0), 16.0); \
  EXPECT_DOUBLE_EQ( myR(1), 28.0); \
  EXPECT_DOUBLE_EQ( myR(2), 6.0);  \
  /* also cover non-zero beta */ \
  constexpr auto beta1 = static_cast<double>(1); \
  pressio::ops::product(pressio::nontranspose(), alpha1, M, VECIN, beta1, myR); \
  EXPECT_DOUBLE_EQ( myR(0), 32.0); \
  EXPECT_DOUBLE_EQ( myR(1), 56.0); \
  EXPECT_DOUBLE_EQ( myR(2), 12.0);  \
  /* alpha = 0, beta = 1 */ \
  constexpr auto alpha0 = static_cast<double>(0);  \
  M(0, 0) = std::nan("0"); /* simulate NaN in input A matrix */ \
  pressio::ops::product(pressio::nontranspose(), alpha0, M, VECIN, beta1, myR); \
  EXPECT_DOUBLE_EQ( myR(0), 32.0); \
  EXPECT_DOUBLE_EQ( myR(1), 56.0); \
  EXPECT_DOUBLE_EQ( myR(2), 12.0); \
  /* alpha = 0, beta = 0 */ \
  myR(0) = std::nan("0"); /* simulate uninitialized NaN */ \
  pressio::ops::product(pressio::nontranspose(), alpha0, M, VECIN, beta0, myR); \
  EXPECT_DOUBLE_EQ( myR(0), 0.0); \
  EXPECT_DOUBLE_EQ( myR(1), 0.0); \
  EXPECT_DOUBLE_EQ( myR(2), 0.0); \

#define OPS_EIGEN_DENSEMATRIX_T_VEC_PROD(VECIN) \
  using M_t = Eigen::MatrixXd; \
  M_t M(4,3);                  \
  M << 1.,0.,2.,2.,1.,3.,0.,0.,1.,2.,3.,4.; \
                               \
  Eigen::VectorXd myR(3);      \
  myR(0) = std::nan("0"); /* simulate uninitialized NaN */ \
  constexpr auto beta0  = static_cast<double>(0); \
  constexpr auto alpha1 = static_cast<double>(1);  \
  pressio::ops::product(pressio::transpose(), alpha1, M, VECIN, beta0, myR); \
  EXPECT_DOUBLE_EQ( myR(0), 14.);  \
  EXPECT_DOUBLE_EQ( myR(1), 11.0); \
  EXPECT_DOUBLE_EQ( myR(2), 8.+6.+6.+12.);  \
  /* also cover non-zero beta */ \
  constexpr auto beta1 = static_cast<double>(1); \
  pressio::ops::product(pressio::transpose(), alpha1, M, VECIN, beta1, myR); \
  EXPECT_DOUBLE_EQ( myR(0), 28.0); \
  EXPECT_DOUBLE_EQ( myR(1), 22.0); \
  EXPECT_DOUBLE_EQ( myR(2), 64.0);  \
  auto y = pressio::ops::product<Eigen::VectorXd>(pressio::transpose(), alpha1, M, VECIN); \
  EXPECT_DOUBLE_EQ( y(0), 14.);  \
  EXPECT_DOUBLE_EQ( y(1), 11.0); \
  EXPECT_DOUBLE_EQ( y(2), 8.+6.+6.+12.);  \
  /* alpha = 0, beta = 1 */ \
  constexpr auto alpha0 = static_cast<double>(0);  \
  M(0, 0) = std::nan("0"); /* simulate NaN in input A matrix */ \
  pressio::ops::product(pressio::transpose(), alpha0, M, VECIN, beta1, myR); \
  EXPECT_DOUBLE_EQ( myR(0), 28.0); \
  EXPECT_DOUBLE_EQ( myR(1), 22.0); \
  EXPECT_DOUBLE_EQ( myR(2), 64.0); \
  /* alpha = 0, beta = 0 */ \
  myR(0) = std::nan("0"); /* simulate uninitialized NaN */ \
  pressio::ops::product(pressio::transpose(), alpha0, M, VECIN, beta0, myR); \
  EXPECT_DOUBLE_EQ( myR(0), 0.0); \
  EXPECT_DOUBLE_EQ( myR(1), 0.0); \
  EXPECT_DOUBLE_EQ( myR(2), 0.0); \


TEST(ops_eigen_level2, dense_matrix_vector_prod)
{
  using V_t = Eigen::VectorXd;
  V_t a(3); a << 4.,2.,6;
  OPS_EIGEN_DENSEMATRIX_VEC_PROD(a);
}

TEST(ops_eigen_level2, dense_matrix_T_vector_prod)
{
  using V_t = Eigen::VectorXd;
  V_t a(4); a << 4.,2.,6,3.;
  OPS_EIGEN_DENSEMATRIX_T_VEC_PROD(a);
}

TEST(ops_eigen_level2, dense_matrix_span_prod)
{
  using V_t = Eigen::VectorXd;
  V_t a(7);
  a(3)=4.;
  a(4)=2.;
  a(5)=6.;

  const auto sp = pressio::span(a, 3, 3);
  OPS_EIGEN_DENSEMATRIX_VEC_PROD(sp);
}

TEST(ops_eigen_level2, dense_matrix_T_span_prod)
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

TEST(ops_eigen_level2, dense_matrix_column_prod){
  using T = Eigen::MatrixXd;
  T a(3,6);
  a.row(0).setConstant(4.);
  a.row(1).setConstant(2.);
  a.row(2).setConstant(6.);
  const auto sp = pressio::column(a, 0);
  OPS_EIGEN_DENSEMATRIX_VEC_PROD(sp);
}

TEST(ops_eigen_level2, dense_matrix_T_column_prod)
{
  using T = Eigen::MatrixXd;
  T a(4,6);
  a.row(0).setConstant(4.);
  a.row(1).setConstant(2.);
  a.row(2).setConstant(6.);
  a.row(3).setConstant(3.);
  const auto sp = pressio::column(a, 0);
  OPS_EIGEN_DENSEMATRIX_T_VEC_PROD(sp);
}

TEST(ops_eigen_level2, dense_matrix_diag_prod)
{
  using T = Eigen::MatrixXd;
  T M0(3,3);
  M0(0,0)=4.;
  M0(1,1)=2.;
  M0(2,2)=6.;

  const auto exp = pressio::diagonal(M0);
  OPS_EIGEN_DENSEMATRIX_VEC_PROD(exp);
}

TEST(ops_eigen_level2, dense_matrix_T_diag_prod)
{
  using T = Eigen::MatrixXd;
  T M0(4,4);
  M0(0,0)=4.;
  M0(1,1)=2.;
  M0(2,2)=6.;
  M0(3,3)=3.;

  const auto exp = pressio::diagonal(M0);
  OPS_EIGEN_DENSEMATRIX_T_VEC_PROD(exp);
}
