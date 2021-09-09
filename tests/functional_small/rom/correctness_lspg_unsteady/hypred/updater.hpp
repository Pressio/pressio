
#ifndef PRESSIO_TEST_ROM_LSPG_UNSTEADY_HYPRED_COMBINER_HPP_
#define PRESSIO_TEST_ROM_LSPG_UNSTEADY_HYPRED_COMBINER_HPP_

#include <gtest/gtest.h>
#include "../../custom_data_types.hpp"

template<class ScalarType>
struct HypRedUpdaterEigen
{
  using vec_operand_type = Eigen::VectorXd;
  using mat_operand_type = Eigen::MatrixXd;
  const std::vector<int> rows_ = {};

  HypRedUpdaterEigen(std::vector<int> rows) : rows_(rows){}

  // a = alpha*a + beta*b (a,b potentially non with same distribution)
  void update_sample_mesh_operand_with_stencil_mesh_one(vec_operand_type & a, ScalarType alpha,
							const vec_operand_type & b, ScalarType beta) const
  {
    for (std::size_t i=0; i<rows_.size(); ++i){
      a(i) = alpha*a(i) + beta*b(rows_[i]);
    }
  }

  void update_sample_mesh_operand_with_stencil_mesh_one(mat_operand_type & a, ScalarType alpha,
							const mat_operand_type & b, ScalarType beta) const
  {
    for (std::size_t i=0; i<rows_.size(); ++i){
      for (int j=0; j<b.cols(); ++j){
	a(i,j) = alpha*a(i,j) + beta*b(rows_[i],j);
      }
    }
  }
};

template<class ScalarType>
struct HypRedUpdaterCustomType
{
  using vec_operand_type = ::pressiotests::MyCustomVector<ScalarType>;
  using mat_operand_type = ::pressiotests::MyCustomMatrix<ScalarType>;
  const std::vector<int> rows_ = {};

  HypRedUpdaterCustomType(std::vector<int> rows) : rows_(rows){}

  // a = alpha*a + beta*b (a,b potentially non with same distribution)
  void update_sample_mesh_operand_with_stencil_mesh_one(vec_operand_type & a, ScalarType alpha,
               const vec_operand_type & b, ScalarType beta) const
  {
    for (std::size_t i=0; i<rows_.size(); ++i){
      a(i) = alpha*a(i) + beta*b(rows_[i]);
    }
  }

  void update_sample_mesh_operand_with_stencil_mesh_one(mat_operand_type & a, ScalarType alpha,
							const mat_operand_type & b, ScalarType beta) const
  {
    for (std::size_t i=0; i<rows_.size(); ++i){
      for (std::size_t j=0; j<b.extent(1); ++j){
	a(i,j) = alpha*a(i,j) + beta*b(rows_[i],j);
      }
    }
  }
};

#endif
