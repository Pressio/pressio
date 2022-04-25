
#include <gtest/gtest.h>
#include "pressio/ops.hpp"

using mat_t = Eigen::MatrixXd;

#define OPS_EIGEN_DENSE_MAT_MAT_PROD(OPERAND) \
  mat_t M(4,3);                     \
  M << 1,0,2, 2,1,3, 0,0,1, 2,2,2;\
  mat_t myR(4,4);         \
  myR(0) = std::nan("0"); /* simulate uninitialized NaN */ \
  constexpr auto beta  = ::pressio::utils::Constants<double>::zero(); \
  constexpr auto alpha = ::pressio::utils::Constants<double>::one();  \
  pressio::ops::product(pressio::nontranspose(), pressio::nontranspose(), \
      alpha, M, OPERAND, beta, myR);\
  EXPECT_DOUBLE_EQ(myR(0,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR(0,1), 3.0); \
  EXPECT_DOUBLE_EQ(myR(0,2), 6.0); \
  EXPECT_DOUBLE_EQ(myR(0,3), 9.0); \
  EXPECT_DOUBLE_EQ(myR(1,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR(1,1), 6.0); \
  EXPECT_DOUBLE_EQ(myR(1,2), 12.0); \
  EXPECT_DOUBLE_EQ(myR(1,3), 18.0); \
  EXPECT_DOUBLE_EQ(myR(2,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR(2,1), 1.0); \
  EXPECT_DOUBLE_EQ(myR(2,2), 2.0); \
  EXPECT_DOUBLE_EQ(myR(2,3), 3.0); \
  EXPECT_DOUBLE_EQ(myR(3,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR(3,1), 6.0); \
  EXPECT_DOUBLE_EQ(myR(3,2), 12.0); \
  EXPECT_DOUBLE_EQ(myR(3,3), 18.0); \
  /* also cover non-zero beta */ \
  constexpr auto beta1 = ::pressio::utils::Constants<double>::one(); \
  pressio::ops::product(pressio::nontranspose(), pressio::nontranspose(), \
      alpha, M, OPERAND, beta1, myR);\
  EXPECT_DOUBLE_EQ(myR(0,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR(0,1), 6.0); \
  EXPECT_DOUBLE_EQ(myR(0,2), 12.0); \
  EXPECT_DOUBLE_EQ(myR(0,3), 18.0); \
  EXPECT_DOUBLE_EQ(myR(1,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR(1,1), 12.0); \
  EXPECT_DOUBLE_EQ(myR(1,2), 24.0); \
  EXPECT_DOUBLE_EQ(myR(1,3), 36.0); \
  EXPECT_DOUBLE_EQ(myR(2,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR(2,1), 2.0); \
  EXPECT_DOUBLE_EQ(myR(2,2), 4.0); \
  EXPECT_DOUBLE_EQ(myR(2,3), 6.0); \
  EXPECT_DOUBLE_EQ(myR(3,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR(3,1), 12.0); \
  EXPECT_DOUBLE_EQ(myR(3,2), 24.0); \
  EXPECT_DOUBLE_EQ(myR(3,3), 36.0); \


#define OPS_EIGEN_DENSE_MAT_T_MAT_PROD(OPERAND) \
  mat_t M(4,3);                     \
  M << 1,0,2, 2,1,3, 0,0,1, 2,2,2;\
                                  \
  mat_t myR(3,3);         \
  myR(0, 0) = std::nan("0"); /* simulate uninitialized NaN */ \
  constexpr auto beta  = ::pressio::utils::Constants<double>::zero(); \
  constexpr auto alpha = ::pressio::utils::Constants<double>::one();  \
  pressio::ops::product(pressio::transpose(), pressio::nontranspose(), \
      alpha, M, OPERAND, beta, myR); \
  EXPECT_DOUBLE_EQ(myR(0,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR(0,1), 5.0); \
  EXPECT_DOUBLE_EQ(myR(0,2), 10.0); \
  EXPECT_DOUBLE_EQ(myR(1,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR(1,1), 3.0); \
  EXPECT_DOUBLE_EQ(myR(1,2), 6.0); \
  EXPECT_DOUBLE_EQ(myR(2,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR(2,1), 8.0); \
  EXPECT_DOUBLE_EQ(myR(2,2), 16.0); \
  /* also cover non-zero beta */ \
  constexpr auto beta1 = ::pressio::utils::Constants<double>::one(); \
  pressio::ops::product(pressio::transpose(), pressio::nontranspose(), \
      alpha, M, OPERAND, beta1, myR); \
  EXPECT_DOUBLE_EQ(myR(0,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR(0,1), 10.0); \
  EXPECT_DOUBLE_EQ(myR(0,2), 20.0); \
  EXPECT_DOUBLE_EQ(myR(1,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR(1,1), 6.0); \
  EXPECT_DOUBLE_EQ(myR(1,2), 12.0); \
  EXPECT_DOUBLE_EQ(myR(2,0), 0.0); \
  EXPECT_DOUBLE_EQ(myR(2,1), 16.0); \
  EXPECT_DOUBLE_EQ(myR(2,2), 32.0); \


#define OPS_EIGEN_DENSE_MAT_T_SELF_PROD \
  mat_t M(4,3);                     \
  M << 1,0,2, 2,1,3, 0,0,1, 2,2,2;\
                                  \
  mat_t myR(3,3);         \
  myR(0, 0) = std::nan("0"); /* simulate uninitialized NaN */ \
  constexpr auto beta  = ::pressio::utils::Constants<double>::zero(); \
  constexpr auto alpha = ::pressio::utils::Constants<double>::one();  \
  pressio::ops::product(pressio::transpose(), pressio::nontranspose(), alpha, M, beta, myR); \
  EXPECT_DOUBLE_EQ(myR(0,0), 9.0); \
  EXPECT_DOUBLE_EQ(myR(0,1), 6.0); \
  EXPECT_DOUBLE_EQ(myR(0,2), 12.0); \
  EXPECT_DOUBLE_EQ(myR(1,0), 6.0); \
  EXPECT_DOUBLE_EQ(myR(1,1), 5.0); \
  EXPECT_DOUBLE_EQ(myR(1,2), 7.0); \
  EXPECT_DOUBLE_EQ(myR(2,0), 12.0); \
  EXPECT_DOUBLE_EQ(myR(2,1), 7.0); \
  EXPECT_DOUBLE_EQ(myR(2,2), 18.0); \
  /* also cover non-zero beta */ \
  constexpr auto beta1 = ::pressio::utils::Constants<double>::one(); \
  pressio::ops::product(pressio::transpose(), pressio::nontranspose(), alpha, M, beta1, myR); \
  EXPECT_DOUBLE_EQ(myR(0,0), 18.0); \
  EXPECT_DOUBLE_EQ(myR(0,1), 12.0); \
  EXPECT_DOUBLE_EQ(myR(0,2), 24.0); \
  EXPECT_DOUBLE_EQ(myR(1,0), 12.0); \
  EXPECT_DOUBLE_EQ(myR(1,1), 10.0); \
  EXPECT_DOUBLE_EQ(myR(1,2), 14.0); \
  EXPECT_DOUBLE_EQ(myR(2,0), 24.0); \
  EXPECT_DOUBLE_EQ(myR(2,1), 14.0); \
  EXPECT_DOUBLE_EQ(myR(2,2), 36.0); \


namespace {
template<class T>
void fillOperand1(T & M)
{
  for (int i=0; i<3; ++i){
    for (int j=0; j<4; ++j){
      M(i,j) = (double) j;
    }
  }
}
template<class T>
void fillOperand2(T & M)
{
  for (int i=0; i<4; ++i){
    for (int j=0; j<3; ++j){
      M(i,j) = (double) j;
    }
  }
}
}//end namespace

TEST(ops_eigen, dense_matrix_dense_matrix_prod)
{
  mat_t operand(3,4);
  fillOperand1(operand);
  OPS_EIGEN_DENSE_MAT_MAT_PROD(operand);
}

TEST(ops_eigen, dense_matrix_T_dense_matrix_prod)
{
  mat_t operand(4,3);
  fillOperand2(operand);
  OPS_EIGEN_DENSE_MAT_T_MAT_PROD(operand);
}

TEST(ops_eigen, dense_matrix_T_self_prod)
{
  OPS_EIGEN_DENSE_MAT_T_SELF_PROD;
}

