
#ifndef PRESSIO_TEST_ROM_LSPG_PRECONDITIONED_DEFAULT_CORRECT_COMMON_PRECOND_HPP_
#define PRESSIO_TEST_ROM_LSPG_PRECONDITIONED_DEFAULT_CORRECT_COMMON_PRECOND_HPP_

#include <gtest/gtest.h>
#include "../custom_data_types.hpp"
#include "pressio/ops.hpp"

struct PreconditionerEigen
{
  using state_type = Eigen::VectorXd;
  using vec_operand_type = Eigen::VectorXd;
  using mat_operand_type = Eigen::MatrixXd;

  PreconditionerEigen(){}

  void operator()(const state_type &, double time, vec_operand_type & operand) const
  {
    std::cout << "EigenPrecond: apply to rank-1 operand\n";
    operand.array() += 1.5;
  }

  void operator()(const state_type &, double time, mat_operand_type & operand) const
  {
    std::cout << "EigenPrecond: apply to rank-2 operand\n";
    operand.array() += 1.5;
  }
};

template<class ScalarType>
struct PreconditionerCustomTypes
{
  using state_type = ::pressiotests::MyCustomVector<ScalarType>;
  using vec_operand_type = ::pressiotests::MyCustomVector<ScalarType>;
  using mat_operand_type = ::pressiotests::MyCustomMatrix<ScalarType>;

  PreconditionerCustomTypes(){}

  void operator()(const state_type &, double time, vec_operand_type & operand) const
  {
    std::cout << "CustomTypesPrecond: apply to rank-1 operand\n";
    for (std::size_t i=0; i<operand.extent(0); ++i){
     operand(i) += 1.5;
    }
  }

  void operator()(const state_type &, double time, mat_operand_type & operand) const
  {
    std::cout << "CustomTypesPrecond: apply to rank-2 operand\n";
    for (std::size_t i=0; i< operand.extent(0); ++i){
      for (std::size_t j=0; j< operand.extent(1); ++j){
        operand(i,j) += 1.5;
      }
    }
  }
};

#endif
