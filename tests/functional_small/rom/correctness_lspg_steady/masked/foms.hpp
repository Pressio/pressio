
#ifndef PRESSIO_TEST_ROM_LSPG_STEADY_MASKED_CORRECT_COMMON_FOM_HPP_
#define PRESSIO_TEST_ROM_LSPG_STEADY_MASKED_CORRECT_COMMON_FOM_HPP_

#include <gtest/gtest.h>
#include "../../custom_data_types.hpp"
#include "pressio/ops.hpp"

struct TrivialFomSteadyEigen
{
  using scalar_type       = double;
  using state_type        = Eigen::VectorXd;
  using residual_type     = state_type;
  int N_ = {};
  const std::vector<int> indices_to_corrupt_ = {};

  TrivialFomSteadyEigen(int N, std::vector<int> ind)
    : N_(N), indices_to_corrupt_(ind){}

  residual_type createResidual() const{ return residual_type(N_); }

  template<class OperandType>
  OperandType createApplyJacobianResult(const OperandType & B) const
  {
    OperandType A(B.rows(), B.cols());
    return A;
  }

  void residual(const state_type & u, residual_type & r) const
  {
    EXPECT_TRUE(u.size()==r.size());
    EXPECT_TRUE(u.size()==N_);

    for (auto i=0; i<r.rows(); ++i){
     r(i) = u(i) + 1.;
    }
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     r(it) = -1114;
    }
  }

  template<class OperandType>
  void applyJacobian(const state_type & state,
                     const OperandType & B,
                     OperandType & A) const
  {
    A = B;
    A.array() += 1.;
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     A.row(it).setConstant(-1114);
    }
  }
};

struct TrivialFomSteadyCustomTypes
{
  using scalar_type       = double;
  using state_type        = ::pressiotests::MyCustomVector<scalar_type>;
  using residual_type     = state_type;
  int N_ = {};
  const std::vector<int> indices_to_corrupt_ = {};

  TrivialFomSteadyCustomTypes(int N, std::vector<int> ind)
    : N_(N), indices_to_corrupt_(ind){}

  residual_type createResidual() const{ return residual_type(N_); }

  template<class OperandType>
  OperandType createApplyJacobianResult(const OperandType & B) const
  {
    OperandType A(N_, B.extent(1));
    return A;
  }

  void residual(const state_type & u, residual_type & r) const
  {
    EXPECT_TRUE((std::size_t)u.extent(0)==(std::size_t)r.extent(0));
    EXPECT_TRUE((std::size_t)u.extent(0)==(std::size_t)N_);

    for (std::size_t i=0; i<r.extent(0); ++i){
     r(i) = u(i) + 1.;
    }
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     r(it) = -1114;
    }
  }

  template<class OperandType>
  void applyJacobian(const state_type & state,
                     const OperandType & B,
                     OperandType & A) const
  {
    A = B;
    for (std::size_t i=0; i< A.extent(0); ++i){
      for (std::size_t j=0; j< A.extent(1); ++j){
        A(i,j) += 1.;
      }
    }
    for (auto & it : indices_to_corrupt_){
      for (std::size_t j=0; j< A.extent(1); ++j){
        A(it,j) = -4232;
      }
    }
  }
};

#endif
