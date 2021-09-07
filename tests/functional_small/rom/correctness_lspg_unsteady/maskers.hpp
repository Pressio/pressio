
#ifndef PRESSIO_TEST_ROM_LSPG_MASKED_STEADY_CORRECT_COMMON_MASKERS_HPP_
#define PRESSIO_TEST_ROM_LSPG_MASKED_STEADY_CORRECT_COMMON_MASKERS_HPP_

#include <gtest/gtest.h>
#include "../custom_data_types.hpp"
#include "pressio/ops.hpp"

// for implicit, masker acts on FOM velicity and FOM apply jac result
struct MaskerSteadyEigen
{
  const std::vector<int> sample_indices_ = {};
  using vec_operand_type = Eigen::VectorXd;
  using mat_operand_type = Eigen::MatrixXd;

  MaskerSteadyEigen(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

  vec_operand_type createApplyMaskResult(const vec_operand_type & operand) const{
    return vec_operand_type(sample_indices_.size());
  }

  mat_operand_type createApplyMaskResult(const mat_operand_type & operand) const{
    return mat_operand_type(sample_indices_.size(), pressio::ops::extent(operand,1));
  }

  void operator()(const vec_operand_type & operand, vec_operand_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      result(i) = operand(sample_indices_[i]);
    }
  }

  void operator()(const mat_operand_type & operand, mat_operand_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      for (int j=0; j< ::pressio::ops::extent(operand,1); ++j){
        result(i,j) = operand(sample_indices_[i],j);
      }
    }
  }
};

struct MaskerSteadyCustomTypes
{
  const std::vector<int> sample_indices_ = {};
  using vec_operand_type = ::pressiotests::MyCustomVector<double>;
  using mat_operand_type = ::pressiotests::MyCustomMatrix<double>;

  MaskerSteadyCustomTypes(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

  vec_operand_type createApplyMaskResult(const vec_operand_type & operand) const{
    return vec_operand_type(sample_indices_.size());
  }

  mat_operand_type createApplyMaskResult(const mat_operand_type & operand) const{
    return mat_operand_type(sample_indices_.size(), operand.extent(1));
  }

  void operator()(const vec_operand_type & operand, vec_operand_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      result(i) = operand(sample_indices_[i]);
    }
  }

  void operator()(const mat_operand_type & operand, mat_operand_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      for (int j=0; j<operand.extent(1); ++j){
        result(i,j) = operand(sample_indices_[i],j);
      }
    }
  }
};


#endif
