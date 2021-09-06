
#ifndef PRESSIO_TEST_ROM_LSPG_STEADY_DEFAULT_CORRECT_FOMS_HPP_
#define PRESSIO_TEST_ROM_LSPG_STEADY_DEFAULT_CORRECT_FOMS_HPP_

#include <gtest/gtest.h>
#include "../../custom_data_types.hpp"

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_Map_decl.hpp>
#endif

struct TrivialFomSteadyEigen
{
  using scalar_type       = double;
  using state_type        = Eigen::VectorXd;
  using residual_type     = state_type;
  int N_ = {};

  TrivialFomSteadyEigen(int N): N_(N){}

  residual_type createResidual() const{ return residual_type(N_); }

  template<class OperandType>
  OperandType createApplyJacobianResult(const OperandType & B) const
  {
    OperandType A(N_, B.cols());
    return A;
  }

  void residual(const state_type & u, residual_type & r) const
  {
    EXPECT_TRUE(u.size()==r.size());
    EXPECT_TRUE(u.size()==N_);

    for (auto i=0; i<r.rows(); ++i){
     r(i) = u(i) + 1.;
    }
  }

  template<class OperandType>
  void applyJacobian(const state_type & state,
                     const OperandType & B,
                     OperandType & A) const
  {
    A = B;
    A.array() += 1.;
  }
};

struct TrivialFomSteadyCustomTypes
{
  using scalar_type       = double;
  using state_type        = ::pressiotests::MyCustomVector<scalar_type>;
  using residual_type     = state_type;
  int N_ = {};

  TrivialFomSteadyCustomTypes(int N): N_(N){}

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
  }
};

#endif
