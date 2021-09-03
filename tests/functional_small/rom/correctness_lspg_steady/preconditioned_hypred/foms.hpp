
#ifndef PRESSIO_TEST_ROM_LSPG_HYPRED_CORRECT_COMMON_FOMS_HPP_
#define PRESSIO_TEST_ROM_LSPG_HYPRED_CORRECT_COMMON_FOMS_HPP_

#include <gtest/gtest.h>
#include "../../custom_data_types.hpp"
#include "pressio/ops.hpp"

const std::vector<int> indices = {0,2,4,6,8,10,12,14};

struct TrivialFomSteadyEigen
{
  using scalar_type       = double;
  using state_type        = Eigen::VectorXd;
  using residual_type     = state_type;
  int N_ = {};

  TrivialFomSteadyEigen(int N) : N_(N){
    EXPECT_TRUE((std::size_t)N==(std::size_t)indices.size());     
  }

  residual_type createResidual() const{ return residual_type(N_); }

  template<class OperandType>
  OperandType createApplyJacobianResult(const OperandType & B) const
  {
    OperandType A(N_, B.cols());
    return A;
  }

  void residual(const state_type & u, residual_type & r) const
  {
    EXPECT_TRUE((std::size_t)u.size()==(std::size_t)15);
    EXPECT_TRUE((std::size_t)u.size()!=(std::size_t)r.size());
    EXPECT_TRUE((std::size_t)r.size()==(std::size_t)N_);

    for (std::size_t i=0; i<indices.size(); ++i){
     r(i) = u(indices[i]);
    }
  }

  template<class OperandType>
  void applyJacobian(const state_type & state,
                     const OperandType & B,
                     OperandType & A) const
  {
    EXPECT_TRUE((std::size_t)B.rows()==(std::size_t)15);
    EXPECT_TRUE((std::size_t)A.rows()==(std::size_t)N_);
    EXPECT_TRUE((std::size_t)state.size()!=(std::size_t)N_);

    for (std::size_t i=0; i<indices.size(); ++i){
      for (int j=0; j< A.cols(); ++j){
        A(i,j) = B(indices[i], j);
      }
    }
  }
};


struct TrivialFomSteadyCustomTypes
{
  using scalar_type       = double;
  using state_type        = ::pressiotests::MyCustomVector<scalar_type>;
  using residual_type     = state_type;
  int N_ = {};

  TrivialFomSteadyCustomTypes(int N) : N_(N){
    EXPECT_TRUE((std::size_t)N==(std::size_t)indices.size());     
  }

  residual_type createResidual() const{ return residual_type(N_); }

  template<class OperandType>
  OperandType createApplyJacobianResult(const OperandType & B) const
  {
    OperandType A(N_, B.extent(1));
    return A;
  }

  void residual(const state_type & u, residual_type & r) const
  {
    EXPECT_TRUE((std::size_t)u.extent(0)==(std::size_t)15);
    EXPECT_TRUE((std::size_t)u.extent(0)!=(std::size_t)r.extent(0));
    EXPECT_TRUE((std::size_t)r.extent(0)==(std::size_t)N_);

    for (std::size_t i=0; i<indices.size(); ++i){
     r(i) = u(indices[i]);
    }
  }

  template<class OperandType>
  void applyJacobian(const state_type & state,
                     const OperandType & B,
                     OperandType & A) const
  {
    EXPECT_TRUE((std::size_t)B.extent(0)==(std::size_t)15);
    EXPECT_TRUE((std::size_t)A.extent(0)==(std::size_t)N_);
    EXPECT_TRUE((std::size_t)state.extent(0)!=(std::size_t)N_);

    for (std::size_t i=0; i<indices.size(); ++i){
      for (int j=0; j< A.extent(1); ++j){
        A(i,j) = B(indices[i], j);
      }
    }
  }
};


#endif
