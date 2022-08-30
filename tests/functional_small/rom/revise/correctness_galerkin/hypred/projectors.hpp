
#ifndef PRESSIO_TEST_ROM_GALERKIN_HYRED_CORRECT_PROJECTORS_COMMON_HPP_
#define PRESSIO_TEST_ROM_GALERKIN_HYRED_CORRECT_PROJECTORS_COMMON_HPP_

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "pressio/ops.hpp"
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_Map_decl.hpp>
#endif

// for explicit, projector acts on FOM velicity
struct ProjectorExplicitEigen
{
  using operator_type = Eigen::MatrixXd; 
  operator_type matrix_;

  ProjectorExplicitEigen(const operator_type & phi) : matrix_(phi){}

  template<class operand_type, class ScalarType, class ResultType>
  void operator()(const operand_type & operand, ScalarType time, ResultType & result) const
  {
    result = matrix_.transpose() * operand;
  }
};

using ProjectorImplicitEigen = ProjectorExplicitEigen;


// for explicit, projector acts on FOM velicity
template<class ScalarType>
struct ProjectorExplicitCustomTypes
{
  using operator_type = ::pressiotests::MyCustomMatrix<ScalarType>;
  operator_type matrix_;

  ProjectorExplicitCustomTypes(const operator_type & phi) : matrix_(phi){}

  // result is the projected RHS, so it is a rom type
  template<class operand_type>
  void operator()(const operand_type & operand, ScalarType time, Eigen::VectorXd & result) const
  {
    // obviously not efficient, just for demonstration
    for (std::size_t k=0; k<matrix_.extent(1); ++k)
    {
      result(k) = 0;
      for (std::size_t i=0; i<matrix_.extent(0); ++i){
        result(k) += matrix_(i,k)*operand(i);
      }
    }
  }
};

// for explicit, projector acts on FOM velicity
template<class ScalarType>
struct ProjectorImplicitCustomTypes
{
  using operator_type = ::pressiotests::MyCustomMatrix<ScalarType>;
  operator_type matrix_;

  ProjectorImplicitCustomTypes(const operator_type & phi) : matrix_(phi){}

  void operator()(const ::pressiotests::MyCustomVector<ScalarType> & operand, 
             ScalarType time, 
             Eigen::VectorXd & result) const
  {
    // obviously not efficient, just for demonstration
    for (std::size_t k=0; k<matrix_.extent(1); ++k)
    {
      result(k) = 0;
      for (std::size_t i=0; i<matrix_.extent(0); ++i){
        result(k) += matrix_(i,k)*operand(i);
      }
    }
  }

  void operator()(const ::pressiotests::MyCustomMatrix<ScalarType> & operand, 
             ScalarType time, 
             Eigen::MatrixXd & result) const
  {
    for (std::size_t i=0; i<matrix_.extent(1); ++i){
      for (std::size_t j=0; j<operand.extent(1); ++j)
      {
        result(i,j) = 0;
        for (std::size_t k=0; k<matrix_.extent(0); ++k){
          result(i,j) += matrix_(k,i)*operand(k,j);
        }
      }
    }
  }
};

#endif
