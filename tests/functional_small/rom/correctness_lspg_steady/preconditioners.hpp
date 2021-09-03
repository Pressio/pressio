
#ifndef PRESSIO_TEST_ROM_LSPG_PRECONDITIONED_DEFAULT_CORRECT_COMMON_PRECOND_HPP_
#define PRESSIO_TEST_ROM_LSPG_PRECONDITIONED_DEFAULT_CORRECT_COMMON_PRECOND_HPP_

#include <gtest/gtest.h>
#include "../custom_data_types.hpp"
#include "pressio/ops.hpp"

struct PreconditionerSteadyEigen
{
  using state_type = Eigen::VectorXd;
  using vec_operand_type = Eigen::VectorXd;
  using mat_operand_type = Eigen::MatrixXd;

  PreconditionerSteadyEigen(){}

  void operator()(const state_type &, vec_operand_type & operand) const
  {
    operand.array() += 1.;
  }

  void operator()(const state_type &, mat_operand_type & operand) const
  {
    operand.array() += 1.;
  }
};

template<class ScalarType>
struct PreconditionerSteadyCustomTypes
{
  using state_type = ::pressiotests::MyCustomVector<ScalarType>;
  using vec_operand_type = ::pressiotests::MyCustomVector<ScalarType>;
  using mat_operand_type = ::pressiotests::MyCustomMatrix<ScalarType>;

  PreconditionerSteadyCustomTypes(){}

  void operator()(const state_type &, vec_operand_type & operand) const
  {
    for (auto i=0; i<operand.extent(0); ++i){
     operand(i) += 1.;
    }
  }

  void operator()(const state_type &, mat_operand_type & operand) const
  {
    for (std::size_t i=0; i< operand.extent(0); ++i){
      for (std::size_t j=0; j< operand.extent(1); ++j){
        operand(i,j) += 1.;
      }
    }
  }
};

#endif
