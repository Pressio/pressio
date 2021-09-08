
#ifndef PRESSIO_TEST_ROM_LSPG_MASKED_UNSTEADY_CORRECT_COMMON_MASKERS_HPP_
#define PRESSIO_TEST_ROM_LSPG_MASKED_UNSTEADY_CORRECT_COMMON_MASKERS_HPP_

#include <gtest/gtest.h>
#include "../custom_data_types.hpp"
#include "pressio/ops.hpp"

struct MaskerEigen
{
  const std::vector<int> sample_indices_ = {};
  using vec_operand_type = Eigen::VectorXd;
  using mat_operand_type = Eigen::MatrixXd;

  MaskerEigen(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

  vec_operand_type createApplyMaskResult(const vec_operand_type & operand) const{
    return vec_operand_type(sample_indices_.size());
  }

  mat_operand_type createApplyMaskResult(const mat_operand_type & operand) const{
    return mat_operand_type(sample_indices_.size(), pressio::ops::extent(operand,1));
  }

  void operator()(const vec_operand_type & operand, double time, vec_operand_type & result) const
  {
    std::cout << "EigenMasker: apply to rank-1 operand\n";
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      result(i) = operand(sample_indices_[i]);
    }
  }

  void operator()(const mat_operand_type & operand, double time, mat_operand_type & result) const
  {
    std::cout << "EigenMasker: apply to rank-2 operand\n";
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      for (int j=0; j< ::pressio::ops::extent(operand,1); ++j){
        result(i,j) = operand(sample_indices_[i],j);
      }
    }
  }
};

struct MaskerCustomTypes
{
  const std::vector<int> sample_indices_ = {};
  using vec_operand_type = ::pressiotests::MyCustomVector<double>;
  using mat_operand_type = ::pressiotests::MyCustomMatrix<double>;

  MaskerCustomTypes(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

  vec_operand_type createApplyMaskResult(const vec_operand_type & operand) const{
    return vec_operand_type(sample_indices_.size());
  }

  mat_operand_type createApplyMaskResult(const mat_operand_type & operand) const{
    return mat_operand_type(sample_indices_.size(), operand.extent(1));
  }

  void operator()(const vec_operand_type & operand, double time, vec_operand_type & result) const
  {
    std::cout << "CustomTypesMasker: apply to rank-1 operand\n";
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      result(i) = operand(sample_indices_[i]);
    }
  }

  void operator()(const mat_operand_type & operand, double time, mat_operand_type & result) const
  {
    std::cout << "CustomTypesMasker: apply to rank-2 operand\n";
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      for (std::size_t j=0; j<operand.extent(1); ++j){
        result(i,j) = operand(sample_indices_[i],j);
      }
    }
  }
};

#endif
